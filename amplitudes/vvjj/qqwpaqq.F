c---------------------------------------------------------------------------
c     lepton helicities have to be summed over externally 
c     ( this routine is for one fixed set of lep. helicities only )
c     
c---------------------------------------------------------------------------
      
      subroutine qqwpaqq(pbar,sign, nlo, L,bos,ps_number,
     1     ans1,ans2,ans3,ans4)
          use globalvars, only: ldoblha
      implicit none
c     
c     Nicolas Kaiser <Nicolas.Kaiser@student.kit.edu>  
c     based on qqwpzqq by Barbara Jaeger
c     Initial version:  July 2013
c     Last modified:     
C     
C     QQWMAQQ calculates the matrix elements**2 for electroweak
c     weak boson pair production by quark quark scattering
C     
C     q1 q3    ---->   q2 q4 W+ A,   W+ ---> f5 f6
C     
c     where the final state particles are 1 charged lepton, 1 neutrino and 1 gamma. 
C     Crossing related processes can be computed as well. Pauli interference 
c     terms for identical fermions are neglected. In pamarticular, only the
c     t-channel exchange of elctroweak bosons is considered. s-channel
c     production of 3 weak bosons is NOT implemented.
c     
C     This code is modified to allow for virtual corrections, more precisely
C     the interference of Born with the finite part of virtual diagrams
C     for 
c     
c     INPUT:  
c     NLO = 0       return uucc = |M_born|^2   etc.
c     for steering the calculation of finite box and pentagon contributions
c     the following additional options are implemented
c     NLO = +4      return cvirt*M_Born contribution
c     remaining cases:   return uucc = 2Re(M_Born^* M_virt)
c     NLO = -4      finite Box contributions only 
c     NLO = -5      true finite Pentagon contributions only
c     
c     integer L (=1,2) is used for the subtraction code to select leptonic 
c     tensors for gluon emission off the upper Line (L=1) or lower quark 
c     line (L=2)
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensorw.inc"
c     
c     electroweak couplings are taken from KOPPLN
c     
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1     V(4,5),A(4,5)
c     alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c     
c     variables for the main part of the program
c     
      double precision  pbar(0:3,4+max_v)
      double precision ans1, ans2, ans3, ans4,res(4),resv(4),
     1     uucs, ddcs, udcc, udss, mgres(4)
      double precision tree(6)
      double precision  p(0:3,4+max_v), p21(0:4), p43(0:4), pew(0:4),
     1     pa(0:4), pwa(0:4)
      integer  sign(4+max_v), nlo, mu, i, j, jj, k, kk,
     1     isig, isig1, isig3,bos,ps_number
      integer  ifl(4,6), js1, js3, L, Lold, is1, is3
      double complex prop21(4), prop43(4)
      double complex mat(4,-1:1,-1:1,9), matv(4,-1:1,-1:1,9,0:2)
      double complex mm(4,-1:1,-1:1),
     1     mv12(4,-1:1,-1:1,0:2), mv34(4,-1:1,-1:1,0:2)
      double complex maw,mzw,mwa,mwz
      double complex rat
      double complex mw(2)
      double complex mbwa(2,3:4,-1:1,-1:1),mbaw(2,3:4,-1:1,-1:1),
     1     mbvwa(2,3:4,-1:1,-1:1,2,0:2),mbvaw(2,3:4,-1:1,-1:1,2,0:2),
     1     zwa(-1:1,-1:1,3),zaw(-1:1,-1:1,3),
     1     zvwa(2,3:4,-1:1,-1:1,3,2,0:2),zvaw(2,3:4,-1:1,-1:1,3,2,0:2),
     1     mbwv,mbvw,gbwv,gbvw,mbvwv,mbvvw,gbvwv,gbvvw
      double complex m5(3,3,3,-1:1,-1:1),mv5(3,3,3,-1:1,-1:1,0:2)
      double complex 
     1     
     1     m1we,m2we,m3we,m4we,w1we,w2we,w3we,w4we,
     2     m1e,m2e,m3e,m4e,
     3     m1aau,m1azu,z1azu,
     4     mw1aau,mw1azu,
     3     m2aau,m2azu,
     4     mw2aau,mw2azu,
     3     m3aau,m3azu,z3azu,
     4     mw3aau,mw3azu,
     3     m4aau,m4azu,
     4     mw4aau,mw4azu,
     2     m1u,m2u,m3u,m4u,
     2     m1awu,z1awu,m1wau,
     2     m2awu,m2wau,
     2     m3awu,z3awu,m3wau,
     2     m4awu,m4wau

      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2) 
      double complex psiwe(2,-1:1,4),psia(2,-1:1,4),
     1     jew(0:5,-1:1,4),ja(0:5,-1:1,4)

      double complex bkjqq(2,-1:1,-1:1,4)
      double complex jw(0:5,-1:1,-1:1,4), 
     2     psiw(2,-1:1,-1:1,4)

      double complex epsww(0:3),
     1     epsaw(0:3),epszw(0:3),epswa(0:3),epswz(0:3),
     1     epsawe(0:3,-1:1,3:4,2), epszwe(0:3,-1:1,3:4,2),
     1     epswae(0:3,-1:1,3:4,2), epswze(0:3,-1:1,3:4,2),
     1     epswee(0:3,-1:1,3:4,2), epswwa(0:3,-1:1,3:4,2),
     1     epszwekk(0:3,-1:1,3:4,2), epswzekk(0:3,-1:1,3:4,2),
     1     epsweekk(0:3,-1:1,3:4,2), epswwakk(0:3,-1:1,3:4,2),
     1     qepswee(3:4),qepszwe(3:4),qepswze(3:4),
     1     qepswwa(3:4)

      double precision fqwe(0:4,4),fqa(0:4,4)
      double precision fq(0:4,4),qee(0:4),qaa(0:4),dummy(0:4)
      double complex zm2i(2:4)
      double complex jj21we(-1:1),jj43we(-1:1),jj21ae(-1:1),jj43ae(-1:1)

      double complex propt1(-1:1,-1:1,4,2), propt2(-1:1,-1:1,4,2),
     1     propw(2),
     1     prop_ee(3), prop_aa(3)
c     
      double complex psiwg,jwg(2,-1:1,-1:1,3),jvwg(2,-1:1,-1:1,3,0:2),ratio
      integer kb,kp
c     
      integer idu,idl
      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c,dotcr
      real*8 theta,mass2
      external contract_Tjj, dotcc, dotrc, dotqj, s1c,dotcr,theta,mass2
      logical ldebug, ldebugv, ldebugvf,linit,lerror
      data linit /.true./, lerror /.false./, ldebug /.false./
      
      logical lmgdebug,lvirtratio
      parameter (lmgdebug = .false.,lvirtratio=.false.)
c     lvirtratio: box- and penlines compute the 1/eps and 1/eps^2 pole 
c     terms of the virt. contributions
c     these terms are proportional to the born: mvirt(div=i)/mborn = const


      parameter (ldebugv=.false., ldebugvf=.false.)
      save ifl, Lold, zm2i, linit
      save jw,psiw,fq
      double complex  zero
      parameter (zero = (0d0,0d0) )


      logical lbornkkterms,lvirtkkterms
      parameter(lbornkkterms=.true.,lvirtkkterms=.true.)
c     
c     variables for virtual corrections
c     

      double precision c2,c2o4pi,scale,pi2o3, cvirtc,
     1     scaleupperline,scalelowerline
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      parameter (pi2o3=pi**2/3d0)
      double precision lnmusqp1p2, lnmusqp3p4
      
      logical lnlo, lbox, lpent, lpt,
     1     lwardtest(40),lward,lwardBox,lwardPent,lwardPent2
      DOUBLE COMPLEX mboxgauge(2,-1:1,-1:1,2,0:2), 
     1     mbox(5,-1:1,-1:1,2,0:2),mboxgaugeborn(2,-1:1,-1:1,2,0:2), 
     2     mboxborn(-1:1,-1:1,2,0:2),
     3     mpentgauge(3,-1:1,-1:1,6,0:2), 
     3     mpent(8,-1:1,-1:1,6,0:2),mpentgaugeborn(3,-1:1,-1:1,6,0:2), 
     4     mpentborn(-1:1,-1:1,6,0:2),jvirtcalc(0:3,-1:1),
     5     wtowakk(0:3)
      double precision  pvirtcalc(0:3,5)

c     color parameters for scale check
      real*8 tf,CF,eulergamma
      parameter (tf=1d0/2d0)
      parameter (eulergamma=5.772156649015328D-1)
      parameter (CF=4.0d0/3.0d0)
      complex*16 cfc,cnf,K_contcheck
      parameter (cnf=(4d0,0d0))
      parameter (cfc=(CF,0d0))

      integer wardIndex,div,divMax,divStep
      wardIndex = 1
      lward = .true.
      lwardBox = .true.
      lwardPent = .true.
      divMax = 0
      divStep = 2
      cvirtc=pi2o3-7d0


      lnlo = NLO.ne.0           ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4 .or. NLO.eq.5 ! call boxline if T
      lpt  = NLO.eq.1 .or. NLO.le.-5 ! call penline if T
      lpent= lpt                ! include pentagons if T

      if(lvirtratio) then
         divMax = 2
         divStep = 1
         lbox = .true.
         lpent = .true.
         cvirtc = 0d0
      endif
c-------------------------

      if (bos.ne.31) then
         print*, "subroutine qqwpaqq contains",
     $        " only ME for W+ A production"
         print*, "but you called it with bos = ",bos
         stop
      endif	
      
c-------------------------
c     
c     
c     define flavors of external quarks for the 4 subprocesses
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
         print*," entering qqwpaqq "
      endif
      do k = 1,4
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,9
                  mat(k,isig1,isig3,i) = 0
                  do div = 0,divMax,divStep
                     matv(k,isig1,isig3,i,div) = 0
                  enddo
               enddo
            enddo
         enddo
      enddo

c     
c     identify fermion line sign factors
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
            print*," are not implemented yet in qqwpaqq"
            stop
         endif
      endif


c     
c     define the internal momenta
c     
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
	 
         p21(mu) =   p(mu,2) - p(mu,1)
         p43(mu) =   p(mu,4) - p(mu,3)
         pew(mu) =   p(mu,6) - p(mu,5)
         pa(mu)  =   p(mu,7)
         pwa(mu) =   pew(mu) + pa(mu)	 	 
	 
      enddo
      p21(4) = mass2(p21)
      p43(4) = mass2(p43)
      pew(4) = mass2(pew)
      pa(4)  = mass2(pa)
      pwa(4) = mass2(pwa)
      
      if(lbox.or.lpent) then
         scaleupperline = -1d0* p21(4)
         scalelowerline = -1d0* p43(4)
      endif

      if (lwarn) then
         do mu=0,3
            dummy(mu) = p(mu,1)+p(mu,2)-
     1           (p(mu,3)+p(mu,4)+p(mu,5)+p(mu,6)+p(mu,7))
         enddo
         if (abs(dummy(0)).gt.1d-8 .or.
     1        abs(dummy(1)).gt.1d-8 .or.
     2        abs(dummy(2)).gt.1d-8 .or.
     3        abs(dummy(3)).gt.1d-8 ) then
            print*," momentum mismatch in qqwpaqq "
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
c     get the vector boson propagator factors
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
c     for box-box and BV graphs we need the propagators for t-channel bosons between quark lines
c     as seen from upper line these W momenta are INCOMING. They are OUTGOINg as seen from lower 
c     line
      do mu = 0,3
         qee(mu) = pew(mu)+p(mu,2) - p(mu,1) ! e ve  emitted on upper line
         qaa(mu) = pa(mu)+p(mu,2) - p(mu,1) !A emitted on upper line
      enddo
      qee(4) = mass2(qee)
      qaa(4) = mass2(qaa)

      prop_ee(1) = 1d0/qee(4)
      prop_ee(2) = 1d0/dcmplx(qee(4)-xm2(2),xmg(2))
      prop_ee(3) = 1d0/dcmplx(qee(4)-xm2(3),xmg(3))
      prop_aa(1) = 1d0/qaa(4)
      prop_aa(2) = 1d0/dcmplx(qaa(4)-xm2(2),xmg(2))
      prop_aa(3) = 1d0/dcmplx(qaa(4)-xm2(3),xmg(3))
      
c     
c     get the external quark spinors (including factor sqrt(2E) )
c     
      call psi0m(4,pbar(0,1),sign(1),psi)
c     
c     get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c     
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2))
      
c     
c     -------------------------------------------------------------------
c     
c     contract with vvtowa tensors to get Vertex-Vertex scattering diagrams
c     
      do isig = -1,1,2
         maw = contract_Tjj(awtowa(0,0,L),jqq(0,isig,1),jqq(0,-1,2))
         mzw = contract_Tjj(zwtowa(0,0,L),jqq(0,isig,1),jqq(0,-1,2))
         mwa = contract_Tjj(watowa(0,0,L),jqq(0,isig,2),jqq(0,-1,1))
         mwz = contract_Tjj(wztowa(0,0,L),jqq(0,isig,2),jqq(0,-1,1))

         do k = 1,2
            mat(k,isig,-1,1) = 
     1           maw*clr(ifl(1,k),1,isig)*clr(ifl(3,k),3,-1)
     2           + mzw*clr(ifl(1,k),2,isig)*clr(ifl(3,k),3,-1)
            do div=0,divMax,divStep
               matv(k,isig,-1,1,div) = (0d0,0d0)
            enddo
         enddo
         do k = 3,4
            mat(k,-1,isig,1) = 
     3           mwa*clr(ifl(1,k),3,-1)*clr(ifl(3,k),1,isig)
     4           + mwz*clr(ifl(1,k),3,-1)*clr(ifl(3,k),2,isig)
            do div=0,divMax,divStep
               matv(k,-1,isig,1,div) = (0d0,0d0)
            enddo
         enddo
      enddo                     !isig


c     -------------------------------------------------------------------
c     -------------------------------------------------------------------
c     
c     prepare box diagrams: attach A/Z to external spinors
c     
c     isig = +-1   : left- and righthanded spinors coupling to A/Z"s
c     
c     Notation for virtual 2-component spinors and momenta
c     
c     W+ -> e+ve   attached to quark number i: psiwe(*,isig,i) with momentum fqwe(mu,i)
c     A  attached to quark number i: psia(*,isig,i) with momentum fqa(mu,i)
c     
c     the fermion current corresponding to a quark line with the emitted A or W->l nu
c     attached next to quark number i is stored in jlv(mu,isig,i) etc. 
c     
c     For the virtual amlitudes the notation, e.g.     jvlv(mu,isig,is,i)
c     is used for the boxline correction to a quark line with one V attached next to
c     quark #i and a free Lorentz index mu for the second attached EW boson. is=+-1
c     refers to the sign factor of this quark (vs.antiquark line). They are recalculated 
c     only if this quark line sign has not been calculated yet for this phase space point
c     (i.e. lws = .true.)  Otherwise they are taken from saved previous calculation
      
      do i = 1,3,2

c     W couples only to left-handed fermions:
	 isig = -1
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,we,
     1        psiwe(1,isig,i),fqwe(0,i))
	 call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,we,
     1        psiwe(1,isig,i+1),fqwe(0,i+1))
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1        psiwe(1,-1,i),fqwe(0,i), jew(0,-1,i)   )
         call curr6(-1,psiwe(1,-1,i+1),fqwe(0,i+1),
     1        psi(1,-1,i),p(0,i), jew(0,-1,i+1) )

c     A also couples to right-handed fermions:
         do isig = -1,1,2
            
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qa,ea,
     1           psia(1,isig,i),fqa(0,i))

            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qa,ea,
     1           psia(1,isig,i+1),fqa(0,i+1))
         enddo                  !isig
         
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1        psia(1,-1,i),fqa(0,i), ja(0,-1,i)   )

         call curr6(1,psia(1,-1,i+1),fqa(0,i+1),
     1        psi(1,-1,i),p(0,i), ja(0,-1,i+1) )
      enddo                     !i


c     
c     -------------------------------------------------------------------
c     -------------------------------------------------------------------
c     
c     now calculate the Vertex-box diagrams; 
      do k = 1,2
         kk = 3-k
c     k=1 and kk=2 is for "box correction" to upper line
c     k=2 and kk=1 is for "box correction" to lower line


         do isig = -1,1,2
            
!     VW-> e nu: (V = A or Z)
            
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
            
!     WW -> A:
            call contract_T1j(wwtoa(0,0,k,L),jqq(0,isig,kk),epsww) !for WW 
            
            do mu = 0,3         ! graphs with uncontracted W from WW -> A
!     (epsww or from VV->e nu (epsvwe/epswve/epswee) 
!     always: k=3/4 for quark type
               
               epswwa(mu,isig,3,k) = epsww(mu)*clr(3,3,isig) !k=3/4 for quark type (0.0d0,0.0d0)
               epswwa(mu,isig,4,k) = epswwa(mu,isig,3,k)
               
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
     1                 epsaw(mu)*clr(3,1,isig)+epszw(mu)*clr(3,2,isig)  
                  epswee(mu,isig,4,k) =
     1                 epsaw(mu)*clr(4,1,isig)+epszw(mu)*clr(4,2,isig)  
               else             !k=1
                  epswee(mu,isig,3,k) =
     1                 epswa(mu)*clr(3,1,isig)+epswz(mu)*clr(3,2,isig)  
                  epswee(mu,isig,4,k) =
     1                 epswa(mu)*clr(4,1,isig)+epswz(mu)*clr(4,2,isig)  
               endif

            enddo               !mu
         enddo                  !isig
      enddo                     !k
      
      do isig = -1,1,2 	
         
         jj21ae(isig) = dotcc(ea,jqq(0,isig,1))
         jj21we(isig) = dotcc(we,jqq(0,isig,1))

         jj43ae(isig) = dotcc(ea,jqq(0,isig,2))
         jj43we(isig) = dotcc(we,jqq(0,isig,2))

      enddo                     !isig
c     
c     --------------
c     
c     now construct the contribution to the amplitude by current contraction 
c     virtual contributions are assembled in subroutine boxline
!     2 bosons attached to 12 line     
c     

      do isig1 = -1,1,2
         do isig3 = -1,1,2
            
c     gauge term for Z/W propagator:
            do k=3,4            !up/down type
               qepszwe(k) = -dotrc(qaa,epszwe(0,isig3,k,1))*zm2i(2)       
               qepswee(k) = -dotrc(qaa,epswee(0,isig3,k,1))*zm2i(3)     
               qepswwa(k) = -dotrc(qee,epswwa(0,isig3,k,1))*zm2i(3)  
               if(.not.lbornkkterms) then
                  qepszwe(k) = 0d0
                  qepswee(k) = 0d0
                  qepswwa(k) = 0d0
               endif   
            enddo               !k
            
            do k = 1,4          ! e+ ve from upper line 
               
               if (isig1.eq.-1) then
                  if (k.eq.2) then  
                     
                     m1we = dotcc(jew(0,isig1,2),epswwa(0,isig3,ifl(3,k),1))	   
                     w1we = qepswwa(ifl(3,k))*jj21we(isig1)    
                     
                     m1e = (m1we+w1we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
                     m2e = 0d0

                  elseif (k.eq.1) then  
                     
                     m2we = dotcc(jew(0,isig1,1),epswwa(0,isig3,ifl(3,k),1))	   
                     w2we = qepswwa(ifl(3,k))*jj21we(isig1)    
                     
                     m2e = (m2we-w2we)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1)
                     m1e = 0d0

                  else          ! k=3,4 			     
                     
                     m1e = 0d0
                     m2e = 0d0
                     
                  endif         !k
                  
               else
                  m1e  =  0d0
                  m2e  =  0d0
               endif            !isig1
               
               
!     A  from upper line:	   
               if (k.lt.3.and.isig3.eq.-1) then !k=1,2
                  m1aau = dotcc(ja(0,isig1,2),epsawe(0,isig3,ifl(3,k),1))
                  m1azu = dotcc(ja(0,isig1,2),epszwe(0,isig3,ifl(3,k),1))

                  z1azu = qepszwe(ifl(3,k))*jj21ae(isig1)  
                  
                  mw1aau =  m1aau*clr(ifl(1,k),1,isig1)*clr(ifl(2,k),1,isig1)
                  mw1azu = (m1azu+z1azu)*clr(ifl(1,k),2,isig1)*clr(ifl(2,k),1,isig1)
                  m1u = mw1aau+mw1azu
                  
                  m2aau = dotcc(ja(0,isig1,1),epsawe(0,isig3,ifl(3,k),1))
                  m2azu = dotcc(ja(0,isig1,1),epszwe(0,isig3,ifl(3,k),1))
                  
                  mw2aau =  m2aau*clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1)
                  mw2azu = (m2azu-z1azu)*clr(ifl(2,k),2,isig1)*clr(ifl(1,k),1,isig1)
                  
                  m2u = mw2aau+mw2azu                   

               elseif (k.gt.2) then !k=3,4
                  m1awu = dotcc(ja(0,isig1,2),epswee(0,isig3,ifl(3,k),1))
                  z1awu = qepswee(ifl(3,k))*jj21ae(isig1)  
                  
                  m1wau = (m1awu+z1awu)*clr(ifl(1,k),3,isig1)*clr(ifl(2,k),1,isig1)
                  m1u = m1wau 
                  
                  m2awu = dotcc(ja(0,isig1,1),epswee(0,isig3,ifl(3,k),1))
                  
                  m2wau = (m2awu-z1awu)*clr(ifl(2,k),3,isig1)*clr(ifl(1,k),1,isig1)
                  
                  m2u = m2wau                  

               else
                  m1u = 0d0
                  m2u = 0d0	  
               endif            !k
               
               mat(k,isig1,isig3,2) =  m1e-m1u+m2e-m2u
               
            enddo               !k
            
         enddo                  !isig3
      enddo                     ! isig1  

      
#if WITH_NLO
      if(lbox) then
!     W+ from upper line:
         do i = 0,3
            pvirtcalc(i,1) = p(i,1)
            pvirtcalc(i,2) =  -1d0 * p(i,2)
            pvirtcalc(i,3) = -1d0 * pew(i)
            pvirtcalc(i,4) =  qee(i)
         enddo
         scale = scaleupperline

         do k=1,4               
            if (k.lt.3) then    ! k=1,2

               do mu=0,3
                  do isig=-1,1,2
                     epswwakk(mu,isig,ifl(3,k),1)=epswwa(mu,isig,ifl(3,k),1)
     $                    - dotcr(epswwa(0,isig,ifl(3,k),1),qee)
     $                    *qee(mu)*zm2i(3)
                     if(.not.lvirtkkterms) then
                        epswwakk(mu,isig,ifl(3,k),1)=epswwa(mu,isig,ifl(3,k),1)
                     endif
                  enddo               
               enddo

               call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              we,epswwakk(0,-1,ifl(3,k),1),-1,1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1


               do div = 0,divMax,divStep
                  do isig3=-1,1,2
                     isig1=-1
                     matv(k,isig1,isig3,2,div) = matv(k,isig1,isig3,2,div)+ 
     $                    ( (mbox(1,isig1,isig3,k,div)-cvirtc*mboxborn(isig1,isig3,k,0) )
     $                    *clr(ifl(1,k),3,isig1)*clr(ifl(2,k),4,isig1) )
                  enddo
               enddo               

            else                !k=3,4
            endif
         enddo                  !k

!     A from upper line
         do i = 0,3
            pvirtcalc(i,1) = p(i,1)
            pvirtcalc(i,2) = -1d0 * p(i,2) 
            pvirtcalc(i,3) = -1d0 * pa(i)
            pvirtcalc(i,4) = qaa(i) 
         enddo
         scale = scaleupperline

         do k=1,4            
            if (k.lt.3) then    !k=1,2
               call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              ea,epsawe(0,-1,ifl(3,k),1),1,-1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1

               do div=0,divMax,divStep
                  do isig1 = -1,1,2
                     isig3 = -1
                     do i = 1,2
                        matv(k,isig1,isig3,2,div) = matv(k,isig1,isig3,2,div) 
     $                       - ( (mbox(1,isig1,isig3,i,div)-cvirtc*mboxborn(isig1,isig3,i,0) )
     $                       *clr(ifl(2,k),1,isig1)*clr(ifl(1,k),1,isig1) )
                     enddo
                  enddo
               enddo

               do mu=0,3
                  do isig=-1,1,2
                     epszwekk(mu,isig,ifl(3,k),1)=epszwe(mu,isig,ifl(3,k),1)
     $                    - dotcr(epszwe(0,isig,ifl(3,k),1),qaa)
     $                    *qaa(mu)*zm2i(2)  
                     if(.not.lvirtkkterms) then
                        epszwekk(mu,isig,ifl(3,k),1)=epszwe(mu,isig,ifl(3,k),1)
                     endif
                  enddo               
               enddo

               call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              ea,epszwekk(0,-1,ifl(3,k),1),1,-1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1


               do div = 0,divMax,divStep
                  do isig1 = -1,1,2
                     isig3 = -1
                     do i = 1,2
                        matv(k,isig1,isig3,2,div) = matv(k,isig1,isig3,2,div)
     $                       - ( (mbox(1,isig1,isig3,i,div)-cvirtc*mboxborn(isig1,isig3,i,0) )
     $                       *clr(ifl(3-i,k),2,isig1)*clr(ifl(i,k),1,isig1) )                         
                     enddo
                  enddo
               enddo

            elseif (k.gt.2) then !k=3,4

               do mu=0,3
                  do isig=-1,1,2
                     epsweekk(mu,isig,ifl(3,k),1)=epswee(mu,isig,ifl(3,k),1)
     $                    - dotcr(epswee(0,isig,ifl(3,k),1),qaa)
     $                    *qaa(mu)*zm2i(3)  
                     if(.not.lvirtkkterms) then
                        epsweekk(mu,isig,ifl(3,k),1)=epswee(mu,isig,ifl(3,k),1)
                     endif
                  enddo               
               enddo

               call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              ea,epsweekk(0,-1,ifl(3,k),1),1,1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1

               do div = 0,divMax,divStep
                  do isig1 = -1,1,2
                     do isig3 = -1,1,2
                        do i = 1,2
                           matv(k,isig1,isig3,2,div) = matv(k,isig1,isig3,2,div) 
     $                          - ( (mbox(1,isig1,isig3,i,div)-cvirtc*mboxborn(isig1,isig3,i,0) )
     $                          *clr(ifl(3-i,k),3,isig1)*clr(ifl(i,k),1,isig1))
                        enddo
                     enddo
                  enddo
               enddo


            endif
         enddo                  !k
      endif                     !lbox	
#endif


c------------------------------------------------------------
      
c     repeat the same for 2 bosons attached to 34 line:
c     
c     
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            
c     gauge term for Z/W propagator:
            do k=3,4            !up/down type
               qepswze(k) = dotrc(qee,epswze(0,isig1,k,2))*zm2i(2) 
               qepswee(k) = dotrc(qee,epswee(0,isig1,k,2))*zm2i(3) 
               qepswwa(k) = dotrc(qaa,epswwa(0,isig1,k,2))*zm2i(3)   
               if(.not.lbornkkterms) then
                  qepszwe(k) = 0d0
                  qepswee(k) = 0d0
                  qepswwa(k) = 0d0
               endif 
            enddo               !k
            
            do k = 1,4          ! W from lower line:
               
               if(isig3.eq.-1) then
                  if (k.lt.3) then !k=1,2
                     
                     m3e = 0d0
                     m4e = 0d0 
                     
                  elseif (k.eq.3) then
                     
                     m3e = 0d0
                     m4we = dotcc(jew(0,isig3,3),epswwa(0,isig1,ifl(1,k),2))
                     w4we = qepswwa(ifl(1,k))*jj43we(isig3)   	   	   	   
                     m4e = (m4we-w4we)*clr(ifl(4,k),4,isig3)*clr(ifl(3,k),3,isig3)

                  elseif (k.eq.4) then
                     
                     m4e = 0d0
                     m3we = dotcc(jew(0,isig3,4),epswwa(0,isig1,ifl(1,k),2))
                     w3we = qepswwa(ifl(1,k))*jj43we(isig3)  	   	   	   
                     m3e = (m3we+w3we)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),3,isig3)
                     
                  endif
                  
               else
                  m3e = 0d0
                  m4e = 0d0
               endif            !isig3  
               
               
!     A from lower line:
               
               if (k.lt.3.and.isig3.eq.-1) then	! k =1,2		
                  m3awu = dotcc(ja(0,isig3,4),epswee(0,isig1,ifl(1,k),2))
                  z3awu = qepswee(ifl(1,k))*jj43ae(isig3)  

                  m3wau = (m3awu+z3awu)*clr(ifl(3,k),3,isig3)*clr(ifl(4,k),1,isig3)      
                  m3u = m3wau 
                  
                  m4awu = dotcc(ja(0,isig3,3),epswee(0,isig1,ifl(1,k),2))
                  m4wau = (m4awu-z3awu)*clr(ifl(4,k),3,isig3)*clr(ifl(3,k),1,isig3)
                  m4u = m4wau 
                  
               elseif (k.gt.2) then !k=3,4
                  m3aau = dotcc(ja(0,isig3,4),epswae(0,isig1,ifl(1,k),2))
                  m3azu = dotcc(ja(0,isig3,4),epswze(0,isig1,ifl(1,k),2))

                  z3azu = qepswze(ifl(1,k))*jj43ae(isig3)  
                  
                  mw3aau =  m3aau*clr(ifl(3,k),1,isig3)*clr(ifl(4,k),1,isig3)
                  mw3azu = (m3azu+z3azu)*clr(ifl(3,k),2,isig3)*clr(ifl(4,k),1,isig3)
                  
                  m3u = mw3aau+mw3azu  
                  
                  m4aau = dotcc(ja(0,isig3,3),epswae(0,isig1,ifl(1,k),2))
                  m4azu = dotcc(ja(0,isig3,3),epswze(0,isig1,ifl(1,k),2))
                  
                  mw4aau =  m4aau*clr(ifl(4,k),1,isig3)*clr(ifl(3,k),1,isig3)
                  mw4azu = (m4azu-z3azu)*clr(ifl(4,k),2,isig3)*clr(ifl(3,k),1,isig3)
                  m4u = mw4aau+mw4azu  

               else
                  m3u = 0d0
                  m4u = 0d0
               endif            !k
               
               mat(k,isig1,isig3,3) = m3e-m3u+m4e-m4u
               
            enddo               !k
            
         enddo                  !isig3
      enddo                     ! isig1  


#ifdef WITH_NLO          
      if (lbox) then
         do i = 0,3
            pvirtcalc(i,1) = p(i,3)
            pvirtcalc(i,2) = -1d0 * p(i,4)
            pvirtcalc(i,3) = -1d0 * pew(i)
            pvirtcalc(i,4) = -1d0 * (p(i,3) - p(i,4) - pew(i))
         enddo
         scale = scalelowerline

         do k=1,4               ! W+ from lower line
            if(k.gt.2) then     ! k=3,4
               
               do mu=0,3
                  do isig=-1,1,2
                     epswwakk(mu,isig,ifl(1,k),2)=epswwa(mu,isig,ifl(1,k),2)
     $                    - dotcr(epswwa(0,isig,ifl(1,k),2),pvirtcalc(0,4))
     $                    *pvirtcalc(mu,4)*zm2i(3)
                     if(.not.lvirtkkterms) then
                        epswwakk(mu,isig,ifl(1,k),2)=epswwa(mu,isig,ifl(1,k),2)
                     endif
                  enddo               
               enddo

               call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              we,epswwakk(0,-1,ifl(1,k),2),-1,1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               do div = 0,divMax,divStep
                  do isig1 = -1,1,2
                  isig3 = -1
                  matv(k,isig1,isig3,3,div) =  matv(k,isig1,isig3,3,div) 
     $                 + ( (mbox(1,isig3,isig1,k-2,div)-cvirtc*mboxborn(isig3,isig1,k-2,0))
     $                 *clr(ifl(4,k),7-k,isig3)
     $                 *clr(ifl(3,k),3,isig3) )                  
                  enddo
               enddo
            endif              
         enddo

!     A from lower line
         do i = 0,3
            pvirtcalc(i,1) = p(i,3)
            pvirtcalc(i,2) = -1d0 * p(i,4)
            pvirtcalc(i,3) = -1d0 * pa(i)
            pvirtcalc(i,4) = -1d0 * (p(i,3) - p(i,4) - pa(i))
         enddo
         scale = scalelowerline

         do k=1,4       
            do mu=0,3
               do isig=-1,1,2
                  epsweekk(mu,isig,ifl(1,k),2)=epswee(mu,isig,ifl(1,k),2)
     $                 - dotcr(epswee(0,isig,ifl(1,k),2),pvirtcalc(0,4))
     $                 *pvirtcalc(mu,4)*zm2i(3)  
                  if(.not.lvirtkkterms) then
                     epsweekk(mu,isig,ifl(1,k),2)=epswee(mu,isig,ifl(1,k),2)
                  endif
               enddo               
            enddo
            
            if(k.lt.3) then     !k=1,2
               call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              ea,epsweekk(0,-1,ifl(1,k),2),-1,1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               do div= 0,divMax,divStep
                  do isig1=-1,1,2
                     isig3=-1
                     do i=1,2
                        matv(k,isig1,isig3,3,div) = matv(k,isig1,isig3,3,div) 
     $                       - ((mbox(1,isig3,isig1,i,div)-cvirtc*mboxborn(isig3,isig1,i,0)) 
     $                   *clr(ifl(5-i,k),3,isig3)*clr(ifl(2+i,k),1,isig3) ) 
                     enddo
                  enddo
               enddo

            elseif (k.gt.2) then !k=3,4

               call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              ea,epswae(0,-1,ifl(1,k),2),1,1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               do div = 0,divMax,divStep
                  do isig1=-1,1,2
                     do isig3=-1,1,2
                        do i=1,2
                           matv(k,isig1,isig3,3,div) = matv(k,isig1,isig3,3,div) 
     $                          - ((mbox(1,isig3,isig1,i,div)-cvirtc*mboxborn(isig3,isig1,i,0))
     $                          *clr(ifl(5-i,k),1,isig3)*clr(ifl(2+i,k),1,isig3) )
                        enddo
                     enddo
                  enddo
               enddo

               do mu=0,3
                  do isig=-1,1,2
                     epswzekk(mu,isig,ifl(1,k),2)=epswze(mu,isig,ifl(1,k),2)
     $                    - dotcr(epswze(0,isig,ifl(1,k),2),pvirtcalc(0,4))
     $                    *pvirtcalc(mu,4)*zm2i(2)  
                     if(.not.lvirtkkterms) then
                        epswzekk(mu,isig,ifl(1,k),2)=epswze(mu,isig,ifl(1,k),2)
                     endif
                  enddo               
               enddo

               call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $              pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $              ea,epswzekk(0,-1,ifl(1,k),2),1,1,
     $              scale,divMax,mboxgauge(1,-1,-1,1,0),
     $              mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $              mboxborn(-1,-1,1,0),lwardtest(wardIndex))
               wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               do div = 0,divMax,divStep
                  do isig1=-1,1,2
                     do isig3=-1,1,2
                        do i=1,2
                           matv(k,isig1,isig3,3,div) = matv(k,isig1,isig3,3,div) 
     $                          - ((mbox(1,isig3,isig1,i,div)-cvirtc*mboxborn(isig3,isig1,i,0) )
     $                          *clr(ifl(5-i,k),2,isig3)*clr(ifl(2+i,k),1,isig3) )  
                        enddo
                     enddo
                  enddo
               enddo

            endif
         enddo
      endif                     !lbox
#endif


c     -----------------------------------------------------------------------
c     
c     next come the W-->WA currents attached to the quark lines. For the virtual 
c     corrections the most effective structure is the contraction of two polarization 
c     vectors with one fermion line. First build these effective polarization vectors
c     from the current wtowa(mu)
c     
c     NOTE: the wtowa currents are NOT conserved. Hence it is necessary 
c     to consider q^mu * q^nu/m_W^2  terms in the W boson propagator


      if (lws(js1) .or. L.ne.Lold) then
         i = 1
	 isig = -1              ! only left-handed q-W coupling
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pwa,wtowa,
     1        psiw(1,isig,is1,i),fq(0,i))
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pwa,wtowa,
     1        psiw(1,isig,is1,i+1),fq(0,i+1))
	 
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1        psiw(1,-1,is1,i),fq(0,i), jw(0,-1,is1,i)   )
         call curr6(-1,psiw(1,-1,is1,i+1),fq(0,i+1),
     1        psi(1,-1,i),p(0,i), jw(0,-1,is1,i+1) )
         
      endif
      
      
      if (lws(js3) .or. L.ne.Lold) then
         i = 3
         isig = -1
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pwa,wtowa,
     1        psiw(1,isig,is3,i),fq(0,i))
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pwa,wtowa,
     1        psiw(1,isig,is3,i+1),fq(0,i+1))

         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1        psiw(1,-1,is3,i),fq(0,i), jw(0,-1,is3,i)   )
         call curr6(-1,psiw(1,-1,is3,i+1),fq(0,i+1),
     1        psi(1,-1,i),p(0,i), jw(0,-1,is3,i+1) )
      endif
      
c     extra terms for massive W propagator:      
      psiwg = -dotrc(pwa,wtowa(0))*zm2i(3)*clr(3,3,-1)
c     
c     
      if(.not.lbornkkterms) then
         psiwg = 0d0
      endif 

      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do kb = 1,3         ! boson id      
               jwg(1,isig1,isig3,kb) = psiwg*dotcc(jqq(0,isig1,1),jqq(0,isig3,2))
     1              *prop43(kb)
               jwg(2,isig1,isig3,kb) = psiwg*dotcc(jqq(0,isig1,1),jqq(0,isig3,2))
     1              *prop21(kb)
            enddo	
         enddo
      enddo

      propw(1) = clr(3,3,-1)**2*prop21(3)
      propw(2) = clr(3,3,-1)**2*prop43(3)
      
      do isig3 = -1,1,2
         isig1 = -1		! only left-handed quarks couple to W
         
c---------------------------

c     box correction to upper line: polarization vectors are 
c     jqq(mu,isig3,2)=j43 with momentum    p43        and
c     wtowa(mu)             with momentum    pwa

         mw(1) =      dotcc(jw(0,isig1,is1,2),jqq(0,isig3,2))
         mw(2) =      dotcc(jw(0,isig1,is1,1),jqq(0,isig3,2))

         do k = 3,4
            propt1(isig1,isig3,k,2) = 
     1           clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1) 
     2           + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2)
            propt2(isig1,isig3,k,2) = 
     1           clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop43(1) 
     2           + clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop43(2)
         enddo
         
         
         if (isig3.eq.-1) then  ! only left-handed quarks couple to W
            do k = 1,2
               kk = 3-k
               mat(k,isig1,isig3,4) = propw(2) *
     1              mw(kk)*clr(3,3,isig1) 
     2              - jwg(1,isig1,isig3,3)*clr(3,3,isig1)*
     3              clr(3,3,isig3)*(3-2*k)	
               
            enddo               !k
         endif                  !isig3
         
         
         do k = 3,4 
            mat(k,isig1,isig3,4) = propt2(isig1,isig3,k,2) *
     1           mw(2)*clr(ifl(1,k),3,isig1)   
     1           + propt1(isig1,isig3,k,2) *
     1           mw(1)*clr(ifl(2,k),3,isig1)        
     1           + jwg(1,isig1,isig3,1)*(clr(ifl(1,k),1,isig1)
     1           -clr(ifl(2,k),1,isig1))*clr(ifl(3,k),1,isig3)
     1           + jwg(1,isig1,isig3,2)*(clr(ifl(1,k),2,isig1)
     1           -clr(ifl(2,k),2,isig1))*clr(ifl(3,k),2,isig3)

         enddo
      enddo                     !isig3	    


#ifdef WITH_NLO          
      if(lbox) then
!     box correction to upper line

         do i = 0,3
            pvirtcalc(i,1) = p(i,1)
            pvirtcalc(i,2) = -1d0 * p(i,2)
            pvirtcalc(i,3) = -1d0 * pwa(i)
            pvirtcalc(i,4) = -1d0 * p43(i)
            do isig = -1,1,2
               jvirtcalc(i,isig) = jqq(i,isig,2)
            enddo
         enddo
         scale = scaleupperline

         call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $        pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $        wtowa,jvirtcalc(0,-1),-1,1,
     $        scale,divMax,mboxgauge(1,-1,-1,1,0),
     $        mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $        mboxborn(-1,-1,1,0),lwardtest(wardIndex))
         wardIndex = wardIndex + 1

         do div = 0,divMax,divStep
            isig1 = -1
            isig3 = -1
            do k=1,2
               matv(k,isig1,isig3,4,div) = propw(2) *
     $              ( mbox(1,isig1,isig3,k,div)
     $              -cvirtc*mboxborn(isig1,isig3,k,0) )*clr(3,3,isig1)
            enddo  


            do isig3 = -1,1,2
               isig1 = -1
               do k=3,4
                  matv(k,isig1,isig3,4,div) = propt2(isig1,isig3,k,2)*
     1                 (mbox(1,isig1,isig3,1,div)
     $                 -cvirtc*mboxborn(isig1,isig3,1,0))*clr(ifl(1,k),3,isig1) 
     2                 + propt1(isig1,isig3,k,2) *
     3                 (mbox(1,isig1,isig3,2,div)
     $                 -cvirtc*mboxborn(isig1,isig3,2,0))*clr(ifl(2,k),3,isig1) 
               enddo
            enddo
         enddo                  ! div
         
      endif                     ! lbox
#endif
      
c     ------------------	    
      
      do isig1 = -1,1,2
         isig3 = -1             ! only left-handed quarks couple to W
         
c     box correction to lower line: polarization vectors are 
c     jqq(mu,isig1,1)=j21 with momentum    p21        and
c     aww/zww(mu)         with momentum    pzz

         mw(1) = dotcc(jw(0,isig3,is3,4),jqq(0,isig1,1))
         mw(2) = dotcc(jw(0,isig3,is3,3),jqq(0,isig1,1))


         do k = 1,2
            propt1(isig1,isig3,k,1) = 
     1           clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)*prop21(1) 
     2           + clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)*prop21(2)
            propt2(isig1,isig3,k,1) = 
     1           clr(ifl(1,k),1,isig1)*clr(ifl(4,k),1,isig3)*prop21(1) 
     2           + clr(ifl(1,k),2,isig1)*clr(ifl(4,k),2,isig3)*prop21(2)
         enddo


         do k = 1,2 
            mat(k,isig1,isig3,5) = propt2(isig1,isig3,k,1) *
     1           mw(2)*clr(ifl(3,k),3,isig3)   
     3           + propt1(isig1,isig3,k,1) *
     1           mw(1)*clr(ifl(4,k),3,isig3)   
     1           + jwg(2,isig1,isig3,1)*(clr(ifl(3,k),1,isig3)-clr(ifl(4,k),1,isig3))*
     1           clr(ifl(1,k),1,isig1)
     1           + jwg(2,isig1,isig3,2)*(clr(ifl(3,k),2,isig3)-clr(ifl(4,k),2,isig3))*
     1           clr(ifl(1,k),2,isig1)
         enddo
         

         if (isig1.eq.-1) then  ! only left-handed quarks couple to W
            do k = 3,4
               kk = 5-k
               mat(k,isig1,isig3,5) = propw(1) *
     1              mw(kk)*clr(3,3,isig1)
     1              - jwg(2,isig1,isig3,3)*clr(3,3,isig1)*clr(3,3,isig3)*(7-2*k)
               
            enddo               !k
         endif                  !isig1
      enddo                     !isig1	    


#ifdef WITH_NLO          
      if(lbox) then
!     box correction to lower line
         do i = 0,3
            pvirtcalc(i,1) = p(i,3)
            pvirtcalc(i,2) = -1d0 * p(i,4)
            pvirtcalc(i,3) = -1d0 * pwa(i)
            pvirtcalc(i,4) = -1d0 * p21(i)
            do isig = -1,1,2
               jvirtcalc(i,isig) = jqq(i,isig,1)
            enddo
         enddo
         scale = scalelowerline

         call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $        pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $        wtowa,jvirtcalc(0,-1),1,1,
     $        scale,divMax,mboxgauge(1,-1,-1,1,0),
     $        mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $        mboxborn(-1,-1,1,0),lwardtest(wardIndex))
         wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do div = 0,divMax,divStep
            do isig1 = -1,1,2
               isig3 = -1
               do k = 1,2
                  matv(k,isig1,isig3,5,div) = propt2(isig1,isig3,k,1)*
     1                 ( (mbox(1,isig3,isig1,1,div)-cvirtc*mboxborn(isig3,isig1,1,0) )*clr(ifl(3,k),3,isig3) )
     2                 + propt1(isig1,isig3,k,1) *
     3                 ( (mbox(1,isig3,isig1,2,div)-cvirtc*mboxborn(isig3,isig1,2,0) )*clr(ifl(4,k),3,isig3) )
               enddo
            enddo

            do k = 3,4
               isig3 = -1
               isig1 = -1
               matv(k,isig1,isig3,5,div) = propw(1) *
     1              (mbox(1,isig3,isig1,k-2,div)-cvirtc*mboxborn(isig3,isig1,k-2,0) )
     2              * clr(3,3,isig1)
            enddo  
         enddo

c     kk term calculation

         do i = 0,3
            pvirtcalc(i,1) = p(i,1)
            pvirtcalc(i,2) = -1d0 * p(i,2)
            pvirtcalc(i,3) = -1d0 * pwa(i)
            pvirtcalc(i,4) = -1d0 * p43(i)
            do isig = -1,1,2
               jvirtcalc(i,isig) = jqq(i,isig,2)
            enddo
         enddo
         scale = scaleupperline

         do mu=0,3
            wtowakk(mu) = - dotcr(wtowa(0),pwa(0))*pwa(mu)*zm2i(3)*clr(3,3,-1)
         enddo

         call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $        pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $        wtowakk(0),jvirtcalc(0,-1),1,1,
     $        scale,divMax,mboxgauge(1,-1,-1,1,0),
     $        mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $        mboxborn(-1,-1,1,0),lwardtest(wardIndex))
         wardIndex = wardIndex + 1

         do div=0,divMax,divStep
            do isig1 = -1,1,2
               do isig3 = -1,1,2
                  do kb = 1,3   ! boson id      
                     jvwg(1,isig1,isig3,kb,div) = (mbox(1,isig1,isig3,2,div)
     1                    -cvirtc*mboxborn(isig1,isig3,2,div))
     1                    *prop43(kb)
                     jvwg(2,isig1,isig3,kb,div) = -(mbox(1,isig1,isig3,1,div)
     1                    -cvirtc*mboxborn(isig1,isig3,1,div))
     1                    *prop21(kb)
                     if(.not.lvirtkkterms) then
                        jvwg(1,isig1,isig3,kb,div) = (0d0,0d0)
                        jvwg(2,isig1,isig3,kb,div) = (0d0,0d0)
                     endif
                  enddo	
               enddo
            enddo
         enddo

         do div = 0,divMax,divStep
            isig1 = -1
            isig3 = -1
            do k=1,2
               matv(k,isig1,isig3,4,div) = matv(k,isig1,isig3,4,div)
     $              - jvwg(1,isig1,isig3,3,div)*clr(3,3,isig1)*
     $              clr(3,3,isig3)*(3-2*k)
            enddo  


            do isig3 = -1,1,2
               isig1 = -1
               do k=3,4
                  matv(k,isig1,isig3,4,div) = matv(k,isig1,isig3,4,div)
     1           + jvwg(1,isig1,isig3,1,div)*(clr(ifl(1,k),1,isig1)
     1           -clr(ifl(2,k),1,isig1))*clr(ifl(3,k),1,isig3)
     1           + jvwg(1,isig1,isig3,2,div)*(clr(ifl(1,k),2,isig1)
     1           -clr(ifl(2,k),2,isig1))*clr(ifl(3,k),2,isig3)
               enddo
            enddo
         enddo                  ! div


         do i = 0,3
            pvirtcalc(i,1) = p(i,3)
            pvirtcalc(i,2) = -1d0 * p(i,4)
            pvirtcalc(i,3) = -1d0 * pwa(i)
            pvirtcalc(i,4) = -1d0 * p21(i)
            do isig = -1,1,2
               jvirtcalc(i,isig) = jqq(i,isig,1)
            enddo
         enddo
         scale = scalelowerline

         do mu=0,3
            wtowakk(mu) = - dotcr(wtowa(0),pwa(0))*pwa(mu)*zm2i(3)*clr(3,3,-1)
         enddo

         call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $        pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $        wtowakk(0),jvirtcalc(0,-1),1,1,
     $        scale,divMax,mboxgauge(1,-1,-1,1,0),
     $        mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $        mboxborn(-1,-1,1,0),lwardtest(wardIndex))
         wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!!  or maybe not?!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do div=0,divMax,divStep
            do isig1 = -1,1,2
               do isig3 = -1,1,2
                  do kb = 1,3   ! boson id      
                     jvwg(1,isig1,isig3,kb,div) = (mbox(1,isig1,isig3,2,div)-cvirtc*mboxborn(isig1,isig3,2,div))
     1                    *prop43(kb)
                     jvwg(2,isig1,isig3,kb,div) = -(mbox(1,isig1,isig3,1,div)-cvirtc*mboxborn(isig1,isig3,1,div))
     1                    *prop21(kb)
                     if(.not.lvirtkkterms) then
                        jvwg(1,isig1,isig3,kb,div) = (0d0,0d0)
                        jvwg(2,isig1,isig3,kb,div) = (0d0,0d0)
                     endif
                  enddo	
               enddo
            enddo
         enddo

         do div = 0,divMax,divStep
            do isig1 = -1,1,2
               isig3 = -1
               do k = 1,2
                  matv(k,isig1,isig3,5,div) = matv(k,isig1,isig3,5,div)
     1           + jvwg(2,isig3,isig1,1,div)*(clr(ifl(3,k),1,isig3)-clr(ifl(4,k),1,isig3))*
     1           clr(ifl(1,k),1,isig1)
     1           + jvwg(2,isig3,isig1,2,div)*(clr(ifl(3,k),2,isig3)-clr(ifl(4,k),2,isig3))*
     1           clr(ifl(1,k),2,isig1) 
               enddo
            enddo

            do k = 3,4
               isig3 = -1
               isig1 = -1
               matv(k,isig1,isig3,5,div) = matv(k,isig1,isig3,5,div)
     1              - jvwg(2,isig3,isig1,3,div)*clr(3,3,isig1)*clr(3,3,isig3)*(7-2*k)
            enddo  
         enddo
c     end of kk term calculation
      endif                     !lbox
#endif


c     ------------------------------------------------------------------------------
c     
c     next do the box-box graphs with one boson (W/A) emitted from the 
c     upper and the other from the lower line. 
c     
c     
c     W from upper line, A from lower line:

      do isig3 = -1,1,2
         isig1 = -1             !(only left-handed fermions couple to W)
c     
         do idu = 1,2           ! fermion id for external W coupling
            do idl = 3,4        ! fermion id for external A  coupling
               mbwa(idu,idl,isig1,isig3) = dotcc(jew(0,isig1,idu),
     1              ja(0,isig3,idl))*clr(3,3,isig1)
            enddo               !idl
         enddo                  !idu
      enddo                     !isig3
      
c     W from lower line, A from upper line:
      
      do isig1 = -1,1,2
         isig3 = -1
c     
         do idu = 1,2           ! fermion id for external A coupling
            do idl = 3,4        ! fermion id for external W coupling
               mbaw(idu,idl,isig1,isig3) = dotcc(jew(0,isig3,idl),
     1              ja(0,isig1,idu))*clr(3,3,isig3)
            enddo               !idl
         enddo                  !idu
      enddo                     !isig1
      
c     -------------------------------


c     for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need
      
      i=3                       ! type of exchanged boson (Z,W,A)  (Z/A contributions cancel)
      do isig3 = -1,1,2
         isig1 = -1   
         
         zwa(isig1,isig3,i) = jj21we(isig1)*jj43ae(isig3)*
     1        zm2i(i)*clr(3,3,isig1)

         if(.not.lbornkkterms) then
            zwa(isig1,isig3,i) = 0d0
         endif 
      enddo                     !isig3

      do isig1 = -1,1,2
         isig3 = -1
         
         zaw(isig1,isig3,i) = jj21ae(isig1)*jj43we(isig3)*
     1        zm2i(i)*clr(3,3,isig3)
         
         if(.not.lbornkkterms) then
            zaw(isig1,isig3,i) = 0d0
         endif 
      enddo                     !isig1


#ifdef WITH_NLO          
      if(lbox) then
c     W from upper line, A from lower line, box correction to upper line
         do idl = 3,4
            do i = 0,3
               pvirtcalc(i,1) = p(i,1)
               pvirtcalc(i,2) = -1d0 * p(i,2)
               pvirtcalc(i,3) = -1d0 * pew(i)
               pvirtcalc(i,4) = qee(i)
               do isig = -1,1,2
                  jvirtcalc(i,isig) = ja(i,isig,idl)
               enddo
            enddo
            scale = scaleupperline

            call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           we,jvirtcalc(0,-1),-1,1,
     $           scale,divMax,mboxgauge(1,-1,-1,1,0),
     $           mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $           mboxborn(-1,-1,1,0),lwardtest(wardIndex))
            wardIndex = wardIndex + 1

            do div = 0,divMax,divStep
               do isig3 = -1,1,2
                  isig1 = -1    !(only left-handed fermions couple to W)
                  do idu = 1,2   
                     mbvwa(idu,idl,isig1,isig3,1,div) = (mbox(1,isig1,isig3,idu,div)
     $                    -cvirtc * mboxborn(isig1,isig3,idu,0) )*clr(3,3,isig1)
                  enddo
               enddo 
            enddo 

c     kk term calculation

            do mu = 0,3
               do isig=-1,1,2
                  jvirtcalc(mu,isig) = -dotcr(ja(0,isig,idl),pvirtcalc(0,4))
     $                 *pvirtcalc(mu,4)*zm2i(3)*clr(3,3,isig)
               enddo
            enddo

            call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           we,jvirtcalc(0,-1),-1,-1,
     $           scale,divMax,mboxgauge(1,-1,-1,1,0),
     $           mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $           mboxborn(-1,-1,1,0),lwardtest(wardIndex))
            wardIndex = wardIndex + 1

            j=5-idl
            do div = 0,divMax,divStep
               isig1=-1
               isig3 = -1
               mbox(1,isig1,isig3,j,div) = -mbox(1,isig1,isig3,j,div)
               mboxgauge(1,isig1,isig3,j,div) = -mboxgauge(1,isig1,isig3,j,div)
               mboxborn(isig1,isig3,j,div) = -mboxborn(isig1,isig3,j,div)
               mboxgaugeborn(1,isig1,isig3,j,div) = -mboxgaugeborn(1,isig1,isig3,j,div)
            enddo

            do div = 0,divMax,divStep
               isig3 = -1
               isig1 = -1       !(only left-handed fermions couple to W)
               do idu = 1,2   
                  zvwa(idu,idl,isig1,isig3,3,1,div) = (mbox(1,isig1,isig3,idu,div)
     $                 -cvirtc * mboxborn(isig1,isig3,idu,div) ) !*clr(3,3,isig1)
                  if(.not.lvirtkkterms) then   
                     zvwa(idu,idl,isig1,isig3,3,1,div) = (0d0,0d0)               
                  endif
               enddo
            enddo 
c     end of kk term calculation

         enddo                  ! idl


c     W from upper line, A from lower line, box correction to lower line
         do idu = 1,2
            do i = 0,3
               pvirtcalc(i,1) = p(i,3)
               pvirtcalc(i,2) = -1d0 * p(i,4)
               pvirtcalc(i,3) = -1d0 * pa(i)
               pvirtcalc(i,4) = -1d0 * (p(i,3) - p(i,4) - pa(i) )
               do isig = -1,1,2
                  jvirtcalc(i,isig) = jew(i,isig,idu)
               enddo
            enddo
            scale = scalelowerline

            call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           ea,jvirtcalc(0,-1),1,-1,
     $           scale,divMax,mboxgauge(1,-1,-1,1,0),
     $           mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $           mboxborn(-1,-1,1,0),lwardtest(wardIndex))
            wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            do div = 0,divMax,divStep
               do isig3 = -1,1,2
                  isig1 = -1    !(only left-handed fermions couple to W)
                  do idl = 3,4   
                     mbvwa(idu,idl,isig1,isig3,2,div) = (mbox(1,isig3,isig1,idl-2,div)
     $                    -cvirtc*mboxborn(isig3,isig1,idl-2,0) )*clr(3,3,isig1)
                  enddo
               enddo             
            enddo          

c     kk term calculation           

            do mu = 0,3
               do isig = -1,1,2
                  jvirtcalc(mu,isig) = -dotcr(jew(0,isig,idu),pvirtcalc(0,4))
     $                 *pvirtcalc(mu,4)*zm2i(3)*clr(3,3,isig)
               enddo
            enddo

            call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           ea,jvirtcalc(0,-1),-1,-1,
     $           scale,divMax,mboxgauge(1,-1,-1,1,0),
     $           mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $           mboxborn(-1,-1,1,0),lwardtest(wardIndex))
            wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            j=3-idu
            do div = 0,divMax,divStep
               isig1=-1
               isig3 = -1
               mbox(1,isig1,isig3,j,div) = -mbox(1,isig1,isig3,j,div)
               mboxgauge(1,isig1,isig3,j,div) = -mboxgauge(1,isig1,isig3,j,div)
               mboxborn(isig1,isig3,j,div) = -mboxborn(isig1,isig3,j,div)
               mboxgaugeborn(1,isig1,isig3,j,div) = -mboxgaugeborn(1,isig1,isig3,j,div)
            enddo
            do div = 0,divMax,divStep
               isig3 = -1
               isig1 = -1       !(only left-handed fermions couple to W)
               do idl = 3,4   
                  zvwa(idu,idl,isig1,isig3,3,2,div) = (mbox(1,isig3,isig1,idl-2,div)
     $                 -cvirtc*mboxborn(isig3,isig1,idl-2,div) ) !*clr(3,3,isig1)
                  if(.not.lvirtkkterms) then
                     zvwa(idu,idl,isig1,isig3,3,2,div) = (0d0,0d0)
                  endif
               enddo
            enddo  
c     end of kk term calculation
         enddo                  ! idu

c     A from upper line, W from lower line, box correction to upper line
         do idl = 3,4

            do i = 0,3
               pvirtcalc(i,1) = p(i,1)
               pvirtcalc(i,2) = -1d0 * p(i,2)
               pvirtcalc(i,3) = -1d0 * pa(i)
               pvirtcalc(i,4) = qaa(i)
               do isig = -1,1,2
                  jvirtcalc(i,isig) = jew(i,isig,idl)
               enddo
            enddo
            scale = scaleupperline

            call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           ea,jvirtcalc(0,-1),1,-1,
     $           scale,divMax,mboxgauge,
     $           mbox,mboxgaugeborn,
     $           mboxborn,lwardtest(wardIndex))
            wardIndex = wardIndex + 1

            do div = 0,divMax,divStep
               do isig1 = -1,1,2
                  isig3 = -1
                  do idu = 1,2   
                     mbvaw(idu,idl,isig1,isig3,1,div) = (mbox(1,isig1,isig3,idu,div)
     $                    -cvirtc*mboxborn(isig1,isig3,idu,0) )*clr(3,3,isig3)
                  enddo
               enddo
            enddo

c     kk term calculation
            do mu = 0,3
               do isig = -1,1,2
                  jvirtcalc(mu,isig) = -dotcr(jew(0,isig,idl),pvirtcalc(0,4))
     $                 *pvirtcalc(mu,4)*zm2i(3)*clr(3,3,isig)
               enddo
            enddo

            call BOX_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           ea,jvirtcalc(0,-1),-1,-1,
     $           scale,divMax,mboxgauge,
     $           mbox,mboxgaugeborn,
     $           mboxborn,lwardtest(wardIndex))
            wardIndex = wardIndex + 1

            j=5-idl
            do div = 0,divMax,divStep
               isig1=-1
               isig3 = -1
               mbox(1,isig1,isig3,j,div) = -mbox(1,isig1,isig3,j,div)
               mboxgauge(1,isig1,isig3,j,div) = -mboxgauge(1,isig1,isig3,j,div)
               mboxborn(isig1,isig3,j,div) = -mboxborn(isig1,isig3,j,div)
               mboxgaugeborn(1,isig1,isig3,j,div) = -mboxgaugeborn(1,isig1,isig3,j,div)
            enddo

            do div = 0,divMax,divStep
               isig1 = -1
               isig3 = -1
               do idu = 1,2   
                  zvaw(idu,idl,isig1,isig3,3,1,div) = (mbox(1,isig1,isig3,idu,div)
     $                 -cvirtc*mboxborn(isig1,isig3,idu,div) ) !*clr(3,3,isig3)
                  if(.not.lvirtkkterms) then                  
                     zvaw(idu,idl,isig1,isig3,3,1,div) = (0d0,0d0)
                  endif
               enddo
            enddo
c     end of kk term calculation
         enddo                  ! idl

c     A from upper line, W from lower line, box correction to lower line
         do idu = 1,2
            do i = 0,3
               pvirtcalc(i,1) = p(i,3)
               pvirtcalc(i,2) = -1d0 * p(i,4)
               pvirtcalc(i,3) = -1d0 * pew(i)
               pvirtcalc(i,4) = -1d0 * (p(i,3) - p(i,4) - pew(i) )
               do isig = -1,1,2
                  jvirtcalc(i,isig) = ja(i,isig,idu)
               enddo
            enddo
            scale = scalelowerline

            call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           we,jvirtcalc(0,-1),-1,1,
     $           scale,divMax,mboxgauge(1,-1,-1,1,0),
     $           mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $           mboxborn(-1,-1,1,0),lwardtest(wardIndex))
            wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            do div = 0,divMax,divStep
               do isig1 = -1,1,2
                  isig3 = -1
                  do idl = 3,4  
                     mbvaw(idu,idl,isig1,isig3,2,div) = (mbox(1,isig3,isig1,idl-2,div)
     $                    - cvirtc*mboxborn(isig3,isig1,idl-2,0))*clr(3,3,isig3)
                  enddo
               enddo
            enddo

c     kk term calculation
            do mu = 0,3
               do isig = -1,1,2
                  jvirtcalc(mu,isig) = -dotcr(ja(0,isig,idu),pvirtcalc(0,4))
     $                 *pvirtcalc(mu,4)*zm2i(3)*clr(3,3,isig)
               enddo
            enddo

            call BOX_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $           pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4), 
     $           we,jvirtcalc(0,-1),-1,-1,
     $           scale,divMax,mboxgauge(1,-1,-1,1,0),
     $           mbox(1,-1,-1,1,0),mboxgaugeborn(1,-1,-1,1,0),
     $           mboxborn(-1,-1,1,0),lwardtest(wardIndex))
            wardIndex = wardIndex + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


            j=3-idu
            do div = 0,divMax,divStep
               isig1=-1
               isig3=-1
               mbox(1,isig1,isig3,j,div) = -mbox(1,isig1,isig3,j,div)
               mboxgauge(1,isig1,isig3,j,div) = -mboxgauge(1,isig1,isig3,j,div)
               mboxborn(isig1,isig3,j,div) = -mboxborn(isig1,isig3,j,div)
               mboxgaugeborn(1,isig1,isig3,j,div) = -mboxgaugeborn(1,isig1,isig3,j,div)
            enddo

             do div = 0,divMax,divStep
               isig1 = -1
               isig3 = -1
               do idl = 3,4  
                  zvaw(idu,idl,isig1,isig3,3,2,div) = (mbox(1,isig3,isig1,idl-2,div)
     $                 - cvirtc*mboxborn(isig3,isig1,idl-2,div)) !*clr(3,3,isig3)!*(-1d0)*isig1*isig3
                  if(.not.lvirtkkterms) then
                     zvaw(idu,idl,isig1,isig3,3,2,div) = (0d0,0d0)
                  endif
               enddo
            enddo
c     end of kk term calculation

         enddo                  ! idu
      endif                     !lbox
#endif

c     -----------------------------
c     
c     now construct amplitudes:	
c     
      
      do isig1 = -1,1,2
         isig3 = -1
         
         
         do k = 1,2
            
            if (isig1.eq.-1) then  
               mbwv = (mbwa(k,3,isig1,isig3)*clr(ifl(3,k),1,isig3)+
     1              mbwa(k,4,isig1,isig3)*clr(ifl(4,k),1,isig3)	 
     2              )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3)
               
               gbwv = (zwa(isig1,isig3,3)*
     1              (clr(ifl(3,k),1,isig3)-clr(ifl(4,k),1,isig3))
     2              )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3)
     3              *(3-2*k)    ! +/- for k=1/2	
               
            else
               mbwv = 0d0
               gbwv = 0d0
            endif


            mbvw = 0d0
            gbvw = 0d0
            do kb = 1,2	
               mbvw = mbvw +(
     1              mbaw(1,3,isig1,isig3)*clr(ifl(1,k),1,isig1)*
     2              clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     3              mbaw(2,3,isig1,isig3)*clr(ifl(2,k),1,isig1)*
     4              clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     5              mbaw(1,4,isig1,isig3)*clr(ifl(1,k),1,isig1)*
     6              clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     7              mbaw(2,4,isig1,isig3)*clr(ifl(2,k),1,isig1)*
     8              clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3) !+
     9              )*prop_aa(kb)
               
               gbvw = 0d0
            enddo               !kb
            
            mat(k,isig1,isig3,6) = mbwv+mbvw+gbvw+gbwv
            
         enddo                  !k
      enddo                     !isig1


c     --------------
      
      
      do isig3 = -1,1,2
         isig1 = -1
         do k = 3,4
            
            if (isig3.eq.-1) then  
               mbvw = (mbaw(1,k,isig1,isig3)*clr(ifl(1,k),1,isig1)+
     1              mbaw(2,k,isig1,isig3)*clr(ifl(2,k),1,isig1)  
     2              )*clr(3,3,isig3)*clr(3,3,isig1)*prop_aa(3)
               
               gbvw = (zaw(isig1,isig3,3)*
     1              (clr(ifl(1,k),1,isig1)-clr(ifl(2,k),1,isig1))
     2              )*clr(3,3,isig3)*clr(3,3,isig1)*prop_aa(3)
     3              *(7-2*k)    ! +/- for k=3/4	
               
            else
               mbvw = 0d0
               gbvw = 0d0
            endif

            mbwv = 0d0
            gbwv = 0d0
            do kb = 1,2	
               mbwv = mbwv +(
     1              mbwa(1,3,isig1,isig3)*clr(ifl(3,k),1,isig3)*
     2              clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     3              mbwa(2,3,isig1,isig3)*clr(ifl(4,k),1,isig3)*
     4              clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     5              mbwa(1,4,isig1,isig3)*clr(ifl(3,k),1,isig3)*
     6              clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     7              mbwa(2,4,isig1,isig3)*clr(ifl(4,k),1,isig3)*
     8              clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3)
     9              )*prop_ee(kb)
               
               gbwv = 0d0

            enddo               !kb


            mat(k,isig1,isig3,6) = mbwv+mbvw+gbvw+gbwv
            
         enddo                  !k
      enddo                     !isig
      
      
c     ----------------------
      
#ifdef WITH_NLO          
      if (lbox) then
         
         do div = 0,divMax,divStep
            do isig1 = -1,1,2
               isig3 = -1

               do i = 1,2       ! 1 is for upper line, 2 for lower line QCD correction
                  jj = 3 + 3*i  ! stored in matv(...,6) and matv(...,9) respectively 
                  do k = 1,2
                     
                     if (isig1.eq.-1) then  
                        mbvwv = (mbvwa(k,3,isig1,isig3,i,div)*clr(ifl(3,k),1,isig3)+
     &                       mbvwa(k,4,isig1,isig3,i,div)*clr(ifl(4,k),1,isig3)
     &                       )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3)

                        gbvwv = ( zvwa(k,3,isig1,isig3,3,i,div)*clr(ifl(3,k),1,isig3)
     1                       -zvwa(k,4,isig1,isig3,3,i,div)*clr(ifl(4,k),1,isig3)
     2                       )*clr(3,3,isig1)*clr(3,3,isig3)*prop_ee(3)
     3                       *(3-2*k) ! +/- for k=1/2	


                     else
                        mbvwv = 0d0
                        gbvwv = 0d0
                     endif

                     mbvvw = 0d0
                     gbvvw = 0d0
                     do kb = 1,2	
                        mbvvw = mbvvw +(
     1                       mbvaw(1,3,isig1,isig3,i,div)*clr(ifl(1,k),1,isig1)*
     2                       clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     3                       mbvaw(2,3,isig1,isig3,i,div)*clr(ifl(2,k),1,isig1)*
     4                       clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     5                       mbvaw(1,4,isig1,isig3,i,div)*clr(ifl(1,k),1,isig1)*
     6                       clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     7                       mbvaw(2,4,isig1,isig3,i,div)*clr(ifl(2,k),1,isig1)*
     8                       clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3)
     9                       )*prop_aa(kb)

                     enddo      !kb

                     matv(k,isig1,isig3,jj,div) = mbvwv+mbvvw+gbvvw+gbvwv
                  enddo         !k
               enddo            !i (jj)
            enddo               !isig1


c     --------------	
            
            do isig3 = -1,1,2
               isig1 = -1
               
               do i = 1,2       ! 1 is for upper line, 2 for lower line QCD correction
                  jj = 3 + 3*i  ! stored in matv(...,6) and matv(...,9) respectively 
                  
                  do k = 3,4
                     
                     if (isig3.eq.-1) then  
                        mbvvw = (mbvaw(1,k,isig1,isig3,i,div)*clr(ifl(1,k),1,isig1)+
     &                       mbvaw(2,k,isig1,isig3,i,div)*clr(ifl(2,k),1,isig1)  
     &                       )*clr(3,3,isig3)*clr(3,3,isig1)*prop_aa(3)
                        
                        gbvvw = (
     1                       zvaw(1,k,isig1,isig3,3,i,div)*clr(ifl(1,k),1,isig1)
     1                       -zvaw(2,k,isig1,isig3,3,i,div)*clr(ifl(2,k),1,isig1)
     2                       )*clr(3,3,isig3)*clr(3,3,isig1)*prop_aa(3)
     3                       *(7-2*k) ! +/- for k=3/4	

                     else
                        mbvvw = 0d0
                        gbvvw = 0d0
                     endif

                     mbvwv = 0d0
                     gbvwv = 0d0
                     do kb = 1,2	
                        mbvwv = mbvwv +(
     &                       mbvwa(1,3,isig1,isig3,i,div)*clr(ifl(3,k),1,isig3)*
     &                       clr(ifl(2,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     &                       mbvwa(2,3,isig1,isig3,i,div)*clr(ifl(4,k),1,isig3)*
     &                       clr(ifl(1,k),kb,isig1)*clr(ifl(4,k),kb,isig3)+
     &                       mbvwa(1,4,isig1,isig3,i,div)*clr(ifl(3,k),1,isig3)*
     &                       clr(ifl(2,k),kb,isig1)*clr(ifl(3,k),kb,isig3)+
     &                       mbvwa(2,4,isig1,isig3,i,div)*clr(ifl(4,k),1,isig3)*
     &                       clr(ifl(1,k),kb,isig1)*clr(ifl(3,k),kb,isig3)
     &                       )*prop_ee(kb)

                     enddo      !kb
                     
                     matv(k,isig1,isig3,jj,div) = mbvwv+mbvvw+gbvvw+gbvwv
                  enddo         !k
               enddo            !i (jj)
            enddo               !isig3
         enddo                  !div
      endif                     !lbox
#endif
      
c     
c     ------------------------------------------------------------------------------
c     
c     and now, finally, the pentagon contributions, i.e. three bosons emitted 
c     from the same quark line


c     compute matrix elements for pentagon on upper line:
c     
c     need m5(pos.j43=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 

      do isig3 = -1,1,2 
         isig1 = -1

         call ket2c(psi(1,isig1,1),.true.,p(0,1),isig1,
     1        p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,1),dummy)
         call bra2c(psi(1,isig1,2),.true.,p(0,2),isig1,
     1        p43,jqq(0,isig3,2),bkjqq(1,isig1,isig3,2),dummy)
         

c     eps1=j43 with v1 arbitrary,eps2=v2,eps3=v3
         m5(1,3,1,isig1,isig3) = s1c(psia(1,isig1,2),we,.true.,isig1,bkjqq(1,isig1,isig3,1))

         m5(1,1,3,isig1,isig3) = s1c(psiwe(1,isig1,2),ea,.true.,isig1,bkjqq(1,isig1,isig3,1))

c     eps1=v1,eps2=j43 with v2 arbitrary,eps3=v3
         m5(2,3,1,isig1,isig3) = s1c(psia(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psiwe(1,isig1,1))

         m5(2,1,3,isig1,isig3) = s1c(psiwe(1,isig1,2),jqq(0,isig3,2),.true.,isig1,psia(1,isig1,1))

c     eps1=v1,eps2=v2,eps3=j43 with v3 arbitrary
         m5(3,3,1,isig1,isig3) = s1c(bkjqq(1,isig1,isig3,2),ea,.true.,isig1,psiwe(1,isig1,1))
         
         m5(3,1,3,isig1,isig3) = s1c(bkjqq(1,isig1,isig3,2),we,.true.,isig1,psia(1,isig1,1))

      enddo                     ! isig3

#ifdef WITH_NLO          
      if (lpent) then	
         do i = 0,3
            pvirtcalc(i,1) = -1d0 * p(i,1)
            pvirtcalc(i,2) =  p(i,2) 
            pvirtcalc(i,3) =  pa(i)
            pvirtcalc(i,4) =  pew(i)
            pvirtcalc(i,5) =  p43(i) 
            do isig = -1,1,2
               jvirtcalc(i,isig) = -1d0 * jqq(i,isig,2)
            enddo
         enddo
         scale = scaleupperline

         call PENT_VBF(psi(1,-1,2),psi(1,-1,1),pvirtcalc(0,1),
     $        pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4),
     $        pvirtcalc(0,5),ea,we,jvirtcalc(0,-1),-1,1,
     $        scale,divMax,mpentgauge(1,-1,-1,1,0),mpent(1,-1,-1,1,0),
     $        mpentgaugeborn(1,-1,-1,1,0),mpentborn(-1,-1,1,0),lwardPent)


         do div = 0,divMax,divStep
            do isig3=-1,1,2
               isig1=-1
               mv5(1,3,1,isig1,isig3,div) = mpent(1,isig1,isig3,6,div)-cvirtc*mpentborn(isig1,isig3,6,0)
               mv5(1,1,3,isig1,isig3,div) = mpent(1,isig1,isig3,5,div)-cvirtc*mpentborn(isig1,isig3,5,0)
               mv5(2,3,1,isig1,isig3,div) = mpent(1,isig1,isig3,4,div)-cvirtc*mpentborn(isig1,isig3,4,0)
               mv5(2,1,3,isig1,isig3,div) = mpent(1,isig1,isig3,3,div)-cvirtc*mpentborn(isig1,isig3,3,0)
               mv5(3,3,1,isig1,isig3,div) = mpent(1,isig1,isig3,2,div)-cvirtc*mpentborn(isig1,isig3,2,0)
               mv5(3,1,3,isig1,isig3,div) = mpent(1,isig1,isig3,1,div)-cvirtc*mpentborn(isig1,isig3,1,0)
            enddo
         enddo
      endif                     !lpent
#endif


      do isig3 = -1,1,2 
         isig1 = -1

         kb = 1
         
         if (isig3.eq.-1) then 
            
            k = 1
            
            mat(k,isig1,isig3,7) = mat(k,isig1,isig3,7) + clr(3,3,-1)**3*prop43(3)*(
     1           m5(2,3,kb,isig1,isig3)*clr(ifl(2,k),kb,isig1)
     2           +m5(3,3,kb,isig1,isig3)*clr(7-ifl(1,k),kb,isig1)
     3           +m5(3,kb,3,isig1,isig3)*clr(ifl(1,k),kb,isig1))
            
#ifdef WITH_NLO          
            if (lpent) then 
               do div = 0,divMax,divStep
                  matv(k,isig1,isig3,7,div) = matv(k,isig1,isig3,7,div) + 
     1                 clr(3,3,-1)**3*prop43(3)*(
     2                 mv5(2,3,kb,isig1,isig3,div)*clr(ifl(2,k),kb,isig1)
     3                 +mv5(3,3,kb,isig1,isig3,div)*clr(7-ifl(1,k),kb,isig1)
     4                 +mv5(3,kb,3,isig1,isig3,div)*clr(ifl(1,k),kb,isig1))
               enddo
            endif
#endif

            k = 2
            
            mat(k,isig1,isig3,7) = mat(k,isig1,isig3,7) 
     1           + clr(3,3,-1)**3*prop43(3)*(
     1           m5(1,3,kb,isig1,isig3)*clr(ifl(2,k),kb,isig1)
     2           +m5(1,kb,3,isig1,isig3)*clr(7-ifl(2,k),kb,isig1)
     3           +m5(2,kb,3,isig1,isig3)*clr(ifl(1,k),kb,isig1))
            
#ifdef WITH_NLO          
            if (lpent) then 
               do div = 0,divMax,divStep
                  matv(k,isig1,isig3,7,div) = matv(k,isig1,isig3,7,div) 
     1                 + clr(3,3,-1)**3*prop43(3)*(
     2                 mv5(1,3,kb,isig1,isig3,div)*clr(ifl(2,k),kb,isig1)
     3                 +mv5(1,kb,3,isig1,isig3,div)*clr(7-ifl(2,k),kb,isig1)
     4                 +mv5(2,kb,3,isig1,isig3,div)*clr(ifl(1,k),kb,isig1))
               enddo
            endif               !lpent
#endif
            
         endif                  !isig3
         
         do k= 3,4
            
            do kp = 1,2
               
               mat(k,isig1,isig3,7) = mat(k,isig1,isig3,7) + 
     1              clr(3,3,-1)*prop43(kp)*clr(ifl(3,k),kp,isig3)*(
     2              m5(1,3,kb,isig1,isig3)*clr(ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     3              +m5(1,kb,3,isig1,isig3)*clr(7-ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     4              +m5(2,3,kb,isig1,isig3)*clr(ifl(2,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     5              +m5(2,kb,3,isig1,isig3)*clr(ifl(1,k),kb,isig1)*clr(7-ifl(2,k),kp,isig1)
     6              +m5(3,3,kb,isig1,isig3)*clr(7-ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     7              +m5(3,kb,3,isig1,isig3)*clr(ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1))

#ifdef WITH_NLO          
               if (lpent) then
                  do div = 0,divMax,divStep
                     matv(k,isig1,isig3,7,div) = matv(k,isig1,isig3,7,div) + 
     1                    clr(3,3,-1)*prop43(kp)*clr(ifl(3,k),kp,isig3)*(
     2                    mv5(1,3,kb,isig1,isig3,div)*clr(ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     3                    +mv5(1,kb,3,isig1,isig3,div)*clr(7-ifl(2,k),kb,isig1)*clr(ifl(1,k),kp,isig1)
     4                    +mv5(2,3,kb,isig1,isig3,div)*clr(ifl(2,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     5                    +mv5(2,kb,3,isig1,isig3,div)*clr(ifl(1,k),kb,isig1)*clr(7-ifl(2,k),kp,isig1)
     6                    +mv5(3,3,kb,isig1,isig3,div)*clr(7-ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1)
     7                    +mv5(3,kb,3,isig1,isig3,div)*clr(ifl(1,k),kb,isig1)*clr(ifl(2,k),kp,isig1))
                  enddo
               endif
#endif

            enddo               !kp
         enddo                  !k
      enddo                     !isig3

c     ------------------------------

c     same for lower line:
c     
c     need m5(pos.j21=1:3,vext1=1:3,vext2=1:3), id of internal boson arbitrary 


      do isig1 = -1,1,2 
         isig3 = -1

         call ket2c(psi(1,isig3,3),.true.,p(0,3),isig3,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,isig3,3),dummy)
         call bra2c(psi(1,isig3,4),.true.,p(0,4),isig3,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,isig3,4),dummy)
         

c     eps1=j21 with v1 arbitrary,eps2=v2,eps3=v3
         m5(1,3,1,isig1,isig3) = s1c(psia(1,isig3,4),we,.true.,isig3,bkjqq(1,isig1,isig3,3))
         m5(1,1,3,isig1,isig3) = s1c(psiwe(1,isig3,4),ea,.true.,isig3,bkjqq(1,isig1,isig3,3))
c     eps1=v1,eps2=j21 with v2 arbitrary,eps3=v3
         m5(2,3,1,isig1,isig3) = s1c(psia(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psiwe(1,isig3,3))
         m5(2,1,3,isig1,isig3) = s1c(psiwe(1,isig3,4),jqq(0,isig1,1),.true.,isig3,psia(1,isig3,3))

c     eps1=v1,eps2=v2,eps3=j21 with v3 arbitrary
         m5(3,3,1,isig1,isig3) = s1c(bkjqq(1,isig1,isig3,4),ea,.true.,isig3,psiwe(1,isig3,3))
         m5(3,1,3,isig1,isig3) = s1c(bkjqq(1,isig1,isig3,4),we,.true.,isig3,psia(1,isig3,3))

      enddo                     !isig1

#ifdef WITH_NLO          
      if (lpent) then	
         do i = 0,3
            pvirtcalc(i,1) = p(i,3)
            pvirtcalc(i,2) = -1d0 * p(i,4) 
            pvirtcalc(i,3) = -1d0 * pa(i)
            pvirtcalc(i,4) = -1d0 * pew(i)
            pvirtcalc(i,5) = -1d0 * p21(i) 
            do isig = -1,1,2
               jvirtcalc(i,isig) = -1d0 * jqq(i,isig,1)
            enddo
         enddo
         scale = scalelowerline

         call PENT_VBF(psi(1,-1,4),psi(1,-1,3),pvirtcalc(0,1),
     $        pvirtcalc(0,2),pvirtcalc(0,3),pvirtcalc(0,4),
     $        pvirtcalc(0,5),ea,we,jvirtcalc(0,-1),-1,1,
     $        scale,divMax,mpentgauge(1,-1,-1,1,0),mpent(1,-1,-1,1,0),
     $        mpentgaugeborn(1,-1,-1,1,0),mpentborn(-1,-1,1,0),lwardPent2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Attention! The helicities are switched!!!!! 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do div = 0,divMax,divStep
            do isig1 = -1,1,2 
               isig3 = -1
               mv5(1,3,1,isig1,isig3,div) = mpent(1,isig3,isig1,6,div)-cvirtc*mpentborn(isig3,isig1,6,0)
               mv5(1,1,3,isig1,isig3,div) = mpent(1,isig3,isig1,5,div)-cvirtc*mpentborn(isig3,isig1,5,0)
               mv5(2,3,1,isig1,isig3,div) = mpent(1,isig3,isig1,4,div)-cvirtc*mpentborn(isig3,isig1,4,0)
               mv5(2,1,3,isig1,isig3,div) = mpent(1,isig3,isig1,3,div)-cvirtc*mpentborn(isig3,isig1,3,0)
               mv5(3,3,1,isig1,isig3,div) = mpent(1,isig3,isig1,2,div)-cvirtc*mpentborn(isig3,isig1,2,0)
               mv5(3,1,3,isig1,isig3,div) = mpent(1,isig3,isig1,1,div)-cvirtc*mpentborn(isig3,isig1,1,0)
            enddo               ! isig1
         enddo                  !div
      endif                     !lpent
#endif


      do isig1 = -1,1,2 
         isig3 = -1

         kb = 1

         if (isig1.eq.-1) then 
            
            k = 3
            
            mat(k,isig1,isig3,8) = mat(k,isig1,isig3,8) + clr(3,3,-1)**3*prop21(3)*(
     1           m5(2,3,kb,isig1,isig3)*clr(ifl(4,k),kb,isig3)
     2           +m5(3,3,kb,isig1,isig3)*clr(7-ifl(3,k),kb,isig3)
     3           +m5(3,kb,3,isig1,isig3)*clr(ifl(3,k),kb,isig3))
            
#ifdef WITH_NLO          
            if (lpent) then
               do div = 0,divMax,divStep
                  matv(k,isig1,isig3,8,div) = matv(k,isig1,isig3,8,div) 
     1                 + clr(3,3,-1)**3*prop21(3)*(
     2                 mv5(2,3,kb,isig1,isig3,div)*clr(ifl(4,k),kb,isig3)
     3                 +mv5(3,3,kb,isig1,isig3,div)*clr(7-ifl(3,k),kb,isig3)
     4                 +mv5(3,kb,3,isig1,isig3,div)*clr(ifl(3,k),kb,isig3))
               enddo
            endif
#endif

            k = 4
            
            mat(k,isig1,isig3,8) = mat(k,isig1,isig3,8) + clr(3,3,-1)**3*prop21(3)*(
     1           m5(1,3,kb,isig1,isig3)*clr(ifl(4,k),kb,isig3)
     2           +m5(1,kb,3,isig1,isig3)*clr(7-ifl(4,k),kb,isig3)
     3           +m5(2,kb,3,isig1,isig3)*clr(ifl(3,k),kb,isig3))
            
#ifdef WITH_NLO          
            if (lpent) then
               do div = 0,divMax,divStep
                  matv(k,isig1,isig3,8,div) = matv(k,isig1,isig3,8,div) 
     1                 + clr(3,3,-1)**3*prop21(3)*(
     2                 mv5(1,3,kb,isig1,isig3,div)*clr(ifl(4,k),kb,isig3)
     3                 +mv5(1,kb,3,isig1,isig3,div)*clr(7-ifl(4,k),kb,isig3)
     4                 +mv5(2,kb,3,isig1,isig3,div)*clr(ifl(3,k),kb,isig3))
               enddo
            endif
#endif
         endif                  !isig1
         
         do k= 1,2
            
            do kp = 1,2
               
               mat(k,isig1,isig3,8) = mat(k,isig1,isig3,8) + 
     1              clr(3,3,-1)*prop21(kp)*clr(ifl(1,k),kp,isig1)*(
     2              m5(1,3,kb,isig1,isig3)*clr(ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     3              +m5(1,kb,3,isig1,isig3)*clr(7-ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     4              +m5(2,3,kb,isig1,isig3)*clr(ifl(4,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     5              +m5(2,kb,3,isig1,isig3)*clr(ifl(3,k),kb,isig3)*clr(7-ifl(4,k),kp,isig3)
     6              +m5(3,3,kb,isig1,isig3)*clr(7-ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     7              +m5(3,kb,3,isig1,isig3)*clr(ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3))

               
#ifdef WITH_NLO          
               if (lpent) then 
                  do div = 0,divMax,divStep
                     matv(k,isig1,isig3,8,div) = matv(k,isig1,isig3,8,div) + 
     1                    clr(3,3,-1)*prop21(kp)*clr(ifl(1,k),kp,isig1)*(
     2                    mv5(1,3,kb,isig1,isig3,div)*clr(ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     3                    +mv5(1,kb,3,isig1,isig3,div)*clr(7-ifl(4,k),kb,isig3)*clr(ifl(3,k),kp,isig3)
     4                    +mv5(2,3,kb,isig1,isig3,div)*clr(ifl(4,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     5                    +mv5(2,kb,3,isig1,isig3,div)*clr(ifl(3,k),kb,isig3)*clr(7-ifl(4,k),kp,isig3)
     6                    +mv5(3,3,kb,isig1,isig3,div)*clr(7-ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3)
     7                    +mv5(3,kb,3,isig1,isig3,div)*clr(ifl(3,k),kb,isig3)*clr(ifl(4,k),kp,isig3))
                  enddo
               endif
#endif

            enddo               !kp

         enddo                  !k


      enddo                     !isig1
      

c----------------------------------------------------------------------------------------

#ifdef WITH_NLO          
      if(lbox) then
         do i = 1,wardIndex-1
            lwardBox = lwardBox.and.lwardtest(i)
         enddo                  ! i
         if(.not.lwardBox) then
            boxcountm = boxcountm+1
            boxcountps(ps_number) = boxcountps(ps_number)+1
         endif
         boxcount2m = boxcount2m+1
         boxcount2ps(ps_number) = boxcount2ps(ps_number)+1
      endif

      if(lpent) then
         lwardPent = lwardPent.and.lwardPent2
         if(.not.lwardPent) then
            pentcountm = pentcountm+1
            pentcountps(ps_number) = pentcountps(ps_number)+1
         endif
         pentcount2m = pentcount2m+1
         pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
      endif

      lward = lwardBox .and. lwardPent
#endif


c     ------------------------------------------------------------------------------
c     ------------------------------------------------------------------------------
      
c     sum the graphs, square them and map them onto ansi etc.

c     i = 1	Vertex-Vertex
c     i = 2,3	Box-Vertex
c     i = 4,5	WtoWA
c     i = 6,9	Box-Box
c     i = 7,8	Pentagon

      
      if(lvirtratio.and.lward)then
         do div = 1,2
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
     $                          matv(k,isig1,isig3,i,0).ne.0d0) then
                              print*,"matv != 0, mat = 0"
                              print*,"k,isig1,isig3,i,matv",k,isig1,isig3,i,matv(k,isig1,isig3,i,0)
                           endif                        
                           if(mat(k,isig1,isig3,6).ne.0d0) then
                              print*,"k, i, sig1, sig3, mborn",k,i,isig1,isig3,mat(k,isig1,isig3,6)
                              print*,"mvirt",matv(k,isig1,isig3,i,div)
                              ratio=(matv(k,isig1,isig3,i,div))/(mat(k,isig1,isig3,6)*K_contcheck) 
                              print*,"ratio",ratio
                           endif
                        else
                           if( mat(k,isig1,isig3,i).eq.0d0 .and.
     $                          matv(k,isig1,isig3,i,0).ne.0d0) then
                              print*,"matv != 0, mat = 0"
                              print*,"k,isig1,isig3,i,matv",k,isig1,isig3,i,matv(k,isig1,isig3,i,0)
                           endif
                           if(mat(k,isig1,isig3,i).ne.0d0) then
                              print*,"k, i, sig1, sig3, mborn",k,i,isig1,isig3,mat(k,isig1,isig3,i)
                              print*,"mvirt",matv(k,isig1,isig3,i,div)
                              ratio = matv(k,isig1,isig3,i,div)/(mat(k,isig1,isig3,i)*K_contcheck) 
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


      do k = 1,4
         res(k) = 0
         resv(k) = 0
         tree(k) = 0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,8
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                 mat(k,isig1,isig3,i)
               enddo
               
               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &              + dimag(mm(k,isig1,isig3))**2
               if (lnlo) then
                  do div = 0,divMax,divStep
                     mv12(k,isig1,isig3,div) = 
     1                    +  matv(k,isig1,isig3,2,div) + matv(k,isig1,isig3,4,div) + 
     2                    matv(k,isig1,isig3,6,div) + matv(k,isig1,isig3,7,div) 
                     mv34(k,isig1,isig3,div) =
     1                    + matv(k,isig1,isig3,3,div) + matv(k,isig1,isig3,5,div) +
     2                    matv(k,isig1,isig3,8,div) + matv(k,isig1,isig3,9,div) 
                  enddo


c     add Born type term and multiply by F_q = alphas*C_2/4pi
c     the factor pi^2/3+9/2 for the born term is after adding the subtraction term
c     and the counter term for the renormalization of the pdfs

                  if (nlo.gt.0) then
                     mv12(k,isig1,isig3,0) = als(1,1)*c2o4pi*
     1                    ( mv12(k,isig1,isig3,0) + mm(k,isig1,isig3)*cvirt )
                     mv34(k,isig1,isig3,0) = als(2,1)*c2o4pi*
     1                    ( mv34(k,isig1,isig3,0) + mm(k,isig1,isig3)*cvirt )
                  else
                     mv12(k,isig1,isig3,0) = 
     1                    als(1,1)*c2o4pi*mv12(k,isig1,isig3,0)
                     mv34(k,isig1,isig3,0) = 
     1                    als(2,1)*c2o4pi*mv34(k,isig1,isig3,0)
                  endif

                  if(divMax.eq.2.and.lward) then
                     if( mm(k,isig1,isig3) .ne. 0d0) then
                        rat = (mv12(k,isig1,isig3,2)+mv34(k,isig1,isig3,2))/ 
     $                       (mm(k,isig1,isig3)  + mat(k,isig1,isig3,6) - mat(k,isig1,isig3,1))
                     endif
                  endif

                  resv(k) = resv(k) + 2*dreal(
     1                 mm(k,isig1,isig3)   *
     1                 conjg( mv12(k,isig1,isig3,0)+mv34(k,isig1,isig3,0) )  )
               endif
            enddo
         enddo  

         tree(k) = res(k) * 9D0 
         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*9d0 ! 9 is the color sum factor
         else
            res(k) = resv(k)*9d0 ! 9 is the color sum factor
         endif
      enddo

#ifdef WITH_MADGRAPH
      if (lmgdebug) then
         call mgr_uucs(pbar,sign,mgres(1))
         call mgr_ddcs(pbar,sign,mgres(2))
         call mgr_udcc(pbar,sign,mgres(3))
         call mgr_udss(pbar,sign,mgres(4))

         print*,"mg debug        madgraph",
     $        "             vbfnlo",
     $        "                 vbfnlo/madgraph "
         print*,"uusc",mgres(1),res(1),res(1)*2/9/mgres(1),
     $        abs(abs(res(1)*2/9/mgres(1))-1) !factor 2/9 is for spin/color summation
         print*,"ddsc",mgres(2),res(2),res(2)*2/9/mgres(2),
     $        abs(abs(res(2)*2/9/mgres(2))-1)
         print*,"ducc",mgres(3),res(3),res(3)*2/9/mgres(3),
     $        abs(abs(res(3)*2/9/mgres(3))-1)
         print*,"duss",mgres(4),res(4),res(4)*2/9/mgres(4),
     $        abs(abs(res(4)*2/9/mgres(4))-1)
         print*
         write(*,*)
         read(*,*)
      endif
#endif


      if (lwarn) then
         resv(1) = ans1         !uucs
         resv(2) = ans2         !ddcs
         resv(3) = ans3         !udcc
         resv(4) = ans4         !udss
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
c     eliminate processes with photon virtuality below cutoff
      if ( abs(p21(4)).lt.qsqAmin .or.  
     1     abs(p43(4)).lt.qsqAmin ) then
         do k = 1,4
            res(k) = 0
            tree(k) = 0
         enddo
      endif

c     set matrix elements to 0 for unstable points
      if(lnlo.and..not.lward) then
         res(1) = 0d0          !uucs 
         res(2) = 0d0          !ddcs 
         res(3) = 0d0          !udcc 
         res(4) = 0d0          !udss 
         tree(1) = 0d0          !uucs 
         tree(2) = 0d0          !ddcs 
         tree(3) = 0d0          !udcc 
         tree(4) = 0d0          !udss 
      endif

      ans1 = res(1)             !uucs 
      ans2 = res(2)             !ddcs 
      ans3 = res(3)             !udcc 
      ans4 = res(4)             !udss 
      
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


