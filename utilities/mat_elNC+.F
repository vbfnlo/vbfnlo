********************************************************************************
********************************************************************************
***                                                                          ***
*** mat_elNC+.F                                                              ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This subroutine (and those accompanying it) calculate the squared matrix ***
*** elements for the boxes and pentagons for the neutral currents. Note that ***
*** we're using FormCalc notation.  Much of this is purloined from           ***
*** FormCalc's SquaredME.F.                                                  ***
***                                                                          ***
********************************************************************************
********************************************************************************

      subroutine squaringNC(p,mass,qk,NC,SInv,TInv,SqMatEl,matFAC,debug)

      implicit none


#include "global.inc"

** Particle label and do-loop variable
**        i  =  1         upper incoming
**           =  2         upper outgoing
**           =  3         lower incoming
**           =  4         lower outgoing
**           =  5         higgs outgoing
      integer i


** mass, mass^2, charge, left /right couplings.  Label = particle i
      double precision mass(5), mass2(5), q(5), gL(4), gR(4)
      double precision dgL(4), dgR(4)


** Helicities of the particles
      integer Hel1
      integer Hel(5)
      common /FChelicities/ Hel

** The momentum invariants
      double precision SInv(n_p+n_v,n_p+n_v), TInv(n_p+n_v,n_p+n_v)
      double precision S, T, T14, U, T24, S34


      double precision MatSUN


** Function calculating the loop*tree product
      double precision sumup
      external sumup


** The tree and loop level components, and dummy tree level (with susy factors)
      double complex Ctree(1), Cloop(1)
      double complex CtreeD1(1), CtreeD2(1)

** Individual loop components - for description, see notebook 6, page 113
** NOTE that, by default, the qqV loops (signified with the ending Qxx) are
** not calculated - instead we use an effective qqV coupling as it's faster
      double complex CloopBuP(1), CloopBuN(1), CloopBlP(1), CloopBlN(1)
      double complex CloopP1(1), CloopP2(1)
      double complex CloopQct(1), CloopQlN(1), CloopQlP(1)
      double complex CloopQuN(1), CloopQuP(1)

** The squared matrix element.  The first entry is current, the second is tree
** (1) or loop (2) level 
      double precision SqMatEl(6,2), SubMatEl


** 4-momenta p(mu,i) of the external quarks. 
        double precision p(0:3,5)


** Subtraction term from real corrections
      double precision Dsub


** Neutral current and current parameter.  cur = NC and:
**              NC = 1            u + c -> u + H + c
**                   2            d + s -> d + H + s  
**                   3            d + c -> d + H + c 
**                   4            d + s -> d + H + s   
      integer NC, cur

** qk is the quark/antiquark pattern label. It is dependent on fsign: 
**               qk = 1          quark-quark
**                  = 2          antiquark-quark
**                  = 3          quark-antiquark
**                  = 4          antiquark-antiquark
      integer qk

** Flag combining current and quark pattern
      integer cQK

** switch deciding whether to run qqV matrix elements: false by default
** NOTE that, if you switch this to .true. you should also set qqVmat to true
** in mat_elCC+.F, and comment out the call to qqV_eff in 
** amplitudes/vvjj/qqhqq.F, or you'll be double counting and IR stuff won't
** cancel and it'll all be horrible.
      logical qqVmat
      parameter (qqVmat = .false.)

** switch deciding whether to multiply amplitudes by SUSY factors, and those
** SUSY factors
      logical matFAC
      double complex susyFAC(2)

** Debugging flag
      integer debug


** These declare the integral id's and the parameters
#include "integrals_NC.h"
#include "mssm.inc"
#include "abbr_NC.h"

** Needed for looptools
#include "looptools.h"

      cur = NC
      
      if (debug .gt. 0 .and. NC .eq. 1) then
         write(*,*)'   '
         write(*,*)'NC, qk =', NC, qk
      end if


** Sorting out the SUSY factors, because that's fun.  
      if (matFAC) then
         susyFAC(1) = SBA*Hgs(HiggsType,1) + CBA*Hgs(HiggsType,2)
         if (HiggsType .eq. 1) then
            susyFAC(2) = SBA
         else if (HiggsType .eq. 2) then
            susyFAC(2) = CBA
         end if
      else
         susyFAC(1) = (1D0,0D0)
         susyFAC(2) = (1D0,0D0)
      end if


** The necessary invariants, in FormCalc notation      
      S = SInv(1,2)
      T = TInv(1,3)
      T14 = TInv(1,4)
      U = TInv(2,3)
      T24 = TInv(2,4)
      S34 = SInv(3,4)  


*      if (debug .gt. 0 .and. NC .eq. 1) then
*      if (NC .eq. 1) then
*         write(*,*)'S,T,T14,U,T24,S34 =', S,T,T14,U,T24,S34
*      end if


      MatSUN = 9D0


** Setting the charges and left / right couplings for the specific current-qk 
** combination: see ntbk 7, p3.  Label: particle(i)
      do i = 1, 5
         q(i) = chargeCur(NC,i)
         mass2(i) = mass(i)**2
      end do
      
      if (cur .gt. 4) then
         if (qk .eq. 2 .or. qk .eq. 4) then
            q(2) = chargeCur(cur,1)
            q(1) = chargeCur(cur,2)
         end if
         if (qk .eq. 3 .or. qk .eq. 4) then
            q(3) = chargeCur(cur,4)
            q(4) = chargeCur(cur,3)
         end if
      end if


      do i = 1, 4
         if (q(i) .eq. Qu) then
            gL(i) = gL3
            gR(i) = gR3
            dgL(i) = dgL3
            dgR(i) = dgR3     
         else if (q(i) .eq. Qd) then
            gL(i) = gL4
            gR(i) = gR4
            dgL(i) = dgL4
            dgR(i) = dgR4
         end if
      end do


      if (debug .gt. 1) then
          do i = 1, 4
             write(*,*)'i =', i
             write(*,*)'mass2(i) =', mass2(i)
             write(*,*)'q(i) =', q(i)
             write(*,*)'gL(i) =', gL(i)
             write(*,*)'gR(i) =', gR(i)
          end do
      end if


** The real corrections
      Dsub = (0D0,0D0)
      call realCor(mass,qk,NC,p,q,Dsub,debug)
*      write(*,*)'Dsub =', Dsub


** The integral values, from (abbr1_angle.F) are calculated 

** For the boxes attached to the lower quark line..
      iintL1 = Cget(T,mass2(3),T24,MW2,MW2,mass2(4))
      iintL2 = Cget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -     MH2,T,MW2,MW2,MW2)
      iintL3 = Dget(MH2,T,mass2(3),mass2(3),
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,T24,MW2,MW2,
     -     MW2,mass2(4))
      iintL4 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -     MH2,2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -     mass2(3),T,mass2(3),MW2,MW2,MW2,mass2(4))
      iintL5 = Cget(T24,T,mass2(3),MW2,mass2(4),mass2(4))
      iintL6 = Cget(T24,T,mass2(3),MZ2,mass2(3),mass2(3))
      iintL7 = Dget(MH2,T24,T,
     -     2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -     mass2(3),mass2(3),MW2,MW2,mass2(4),mass2(4))
      iintL8 = Dget(MH2,T24,T,
     -     2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -     mass2(3),mass2(3),MZ2,MZ2,mass2(3),mass2(3))


** For the boxes attached to the upper quark line..
      iintU1 = Cget(MH2,
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,T,MW2,MW2,MW2)
      iintU2 = iintL2
      iintU3 = Dget(MH2,
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(1),
     -     mass2(1),T,S34,MW2,MW2,MW2,mass2(2))
      iintU4 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -     MH2,mass2(1),mass2(1),T,T14,MW2,MW2,MW2,mass2(2))
      iintU5 = Cget(mass2(1),
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,T14,MW2,
     -     mass2(2),mass2(2))
      iintU6 = Cget(mass2(1),
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,T14,MZ2,
     -     mass2(1),mass2(1))
      iintU7 = Dget(MH2,mass2(1),
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(1),S34,
     -     T14,MW2,MW2,mass2(2),mass2(2))
      iintU8 = Dget(MH2,mass2(1),
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(1),S34,
     -     T14,MZ2,MZ2,mass2(1),mass2(1))

      
** For the pentagons ..
      if (cur .eq. 1 .or. cur .eq. 4) then
         if (qk .eq. 1 .or. qk .eq. 4) then
            cQK = 1
         else if (qk .eq. 2 .or. qk .eq. 3) then
            cQK = 2
         end if
      else if (cur .eq. 2 .or. cur .eq. 3) then
         if (qk .eq. 1 .or. qk .eq. 4) then
            cQK = 4
         else if (qk .eq. 2 .or. qk .eq. 3) then
            cQK = 3
         end if
      end if
         
      iintP5 = Eget(MH2,
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(1),
     -     3*mass2(1) + 2*mass2(3) + MH2 - S - T - T14,T24,T,S34,
     -     mass2(3),mass2(1),mass2(3),MZ2,MZ2,MZ2,mass2(1),mass2(3))
      iintP6 = Eget(T,MH2,S34,
     -     3*mass2(1) + 2*mass2(3) + MH2 - S - T - T14,mass2(3),
     -     2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(1),
     -     mass2(3),mass2(1),T24,0.D0,MZ2,MZ2,mass2(1),mass2(3))
      iintP7 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -     MH2,mass2(1),S,mass2(3),T,S34,
     -     2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -     mass2(1),mass2(3),0.D0,MZ2,MZ2,mass2(1),mass2(3))
      iintP8 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -     MH2,mass2(1),S,mass2(3),T,S34,
     -     2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -     mass2(1),mass2(3),MZ2,MZ2,MZ2,mass2(1),mass2(3))
      iintP9 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -     MH2,mass2(1),3*mass2(1) + 2*mass2(3) + MH2 - S34 - T - U,
     -     mass2(3),T,T14,T24,mass2(1),mass2(3),0.D0,MZ2,MZ2,mass2(1),
     -     mass2(3))
      iintP11 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -     U,MH2,mass2(1),U,mass2(3),T,T14,
     -     2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -     mass2(1),mass2(3),0.D0,MZ2,MZ2,mass2(1),mass2(3))
      iintP12 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -     U,MH2,mass2(1),U,mass2(3),T,T14,
     -     2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -     mass2(1),mass2(3),MZ2,MZ2,MZ2,mass2(1),mass2(3))

      if (cQK .eq. 2 .or. cQK .eq. 3) then
         iintP4 = iintP5
         iintP5 = iintP6
         iintP6 = iintP7
         iintP7 = iintP8
         iintP8 = iintP9
         iintP10 = iintP11
         iintP11 = iintP12
      end if

      if (cQK .eq. 1 .or. cQK .eq. 4) then
         iintP1 = Dget(T,MH2,mass2(3),mass2(3),
     -        2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,T24,0.D0,MZ2,
     -        MZ2,mass2(3))
         iintP2 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -        MH2,2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -        mass2(3),T,mass2(3),0.D0,MZ2,MZ2,mass2(3))
         iintP3 = Dget(MH2,
     -        2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(3),T24,
     -        T,mass2(3),MZ2,MZ2,MZ2,mass2(3))
         iintP4 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -        MH2,2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -        mass2(3),T,mass2(3),MZ2,MZ2,MZ2,mass2(3))
         iintP10 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,MH2,mass2(1),3*mass2(1) + 2*mass2(3) + MH2 - S34 - T - 
     -        U,mass2(3),T,T14,T24,mass2(1),mass2(3),MZ2,MZ2,MZ2,
     -        mass2(1),mass2(3))
      end if

      if (cQK .eq. 1 .or. cQK .eq. 3) then
         iintP15 = Eget(MH2,
     -        2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(1),
     -        3*mass2(1) + 2*mass2(3) + MH2 - S - T - T14,T24,T,S34,
     -        mass2(3),mass2(1),mass2(3),MW2,MW2,MW2,mass2(2),mass2(4))
         iintP16 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,MH2,mass2(1),U,mass2(3),T,T14,
     -        2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -        mass2(1),mass2(3),MW2,MW2,MW2,mass2(2),mass2(4))
         if (cQK .eq. 3) then
            iintP12 = iintP15
            iintP13 = iintP16
         end if
      end if

      if (cQK .eq. 1) then
         iintP13 = Dget(MH2,
     -        2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(3),T24,
     -        T,mass2(3),MW2,MW2,MW2,mass2(4))
         iintP14 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,MH2,2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -        mass2(3),T,mass2(3),MW2,MW2,MW2,mass2(4))
      end if

      if (cQK .eq. 2 .or. cQK .eq. 3) then
         iintP1 = Dget(MH2,mass2(1),S,mass2(3),S34,
     -        2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,MZ2,MZ2,
     -        mass2(1),mass2(3))
         iintP2 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -        MH2,T24,mass2(3),T,mass2(3),0.D0,MZ2,MZ2,mass2(3))
         iintP3 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -        MH2,T24,mass2(3),T,mass2(3),MZ2,MZ2,MZ2,mass2(3))
         iintP9 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -        MH2,mass2(1),3*mass2(1) + 2*mass2(3) + MH2 - S34 - T - U,
     -        mass2(3),T,T14,T24,mass2(1),mass2(3),MZ2,MZ2,MZ2,mass2(1),
     -        mass2(3))
      end if

      if (cQK .eq. 2 .or. cQK .eq. 4) then
         iintP14 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,MH2,mass2(1),S,mass2(3),T,S34,
     -        2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,
     -        mass2(1),mass2(3),MW2,MW2,MW2,mass2(2),mass2(4))
         iintP15 = Eget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,MH2,mass2(1),3*mass2(1) + 2*mass2(3) + MH2 - S34 - T - 
     -        U,mass2(3),T,T14,T24,mass2(1),mass2(3),MW2,MW2,MW2,
     -        mass2(2),mass2(4))
         if (cQK .eq. 4) then
            iintP13 = iintP14
            iintP14 = iintP15
         end if
      end if

      if (cQK .eq. 2) then
         iintP12 = Dget(MH2,mass2(1),S,mass2(3),S34,
     -        2*mass2(1) + 2*mass2(3) + 2*MH2 - S34 - T14 - T24,MW2,MW2,
     -        mass2(2),mass2(4))
         iintP13 = Dget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,MH2,T24,mass2(3),T,mass2(3),MW2,MW2,MW2,mass2(4))
      end if
         
      
** For the qqV matrix elements .. (these are called iff no approximations are 
** being considered)
      if (qqVmat) then     
         iintQ1 = Bget(T,mass2(1),mass2(1))
         iintQ2 = Bget(T,mass2(2),mass2(2))
         iintQ3 = Cget(mass2(1),T,mass2(1),0.D0,mass2(1),mass2(1))
         iintQ4 = Cget(mass2(1),T,mass2(1),MW2,mass2(2),mass2(2))
         iintQ5 = Cget(mass2(1),T,mass2(1),MZ2,mass2(1),mass2(1))
         iintQ6 = Bget(T,MW2,MW2)
         iintQ7 = Cget(T,mass2(1),mass2(1),MW2,MW2,mass2(2))
         iintQ8 = Bget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -        mass2(3),mass2(3))
         iintQ9 = Bget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,
     -        mass2(4),mass2(4))
         iintQ10 = Cget(mass2(3),
     -        2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(3),0.D0,
     -        mass2(3),mass2(3))
         iintQ11 = Cget(mass2(3),
     -        2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(3),MW2,
     -        mass2(4),mass2(4))
         iintQ12 = Cget(mass2(3),
     -        2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U,mass2(3),MZ2,
     -        mass2(3),mass2(3))
         iintQ13 = Bget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,MW2,MW2)
         iintQ14 = Cget(2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - 
     -        U,mass2(3),mass2(3),MW2,MW2,mass2(4))
      end if


** Clearing the values of the squared matrix element
      SubMatEl = 0D0
      SqMatEl(NC,1) = 0D0
      SqMatEl(NC,2) = 0D0


** Looping over permissible helicity combinations
      do Hel1 = 1, 4
         if (Hel1 .eq. 1) then
            Hel(1) = -1
            Hel(2) = -1
            Hel(3) = -1
            Hel(4) = -2
            Hel(5) = -1
         else if (Hel1 .eq. 2) then
            Hel(1) = -1
            Hel(2) = 1
            Hel(3) = -1
            Hel(4) = -2
            Hel(5) = 1
         else if (Hel1 .eq. 3) then
            Hel(1) = 1
            Hel(2) = -1
            Hel(3) = 1
            Hel(4) = -2
            Hel(5) = -1
         else if (Hel1 .eq. 4) then
            Hel(1) = 1
            Hel(2) = 1
            Hel(3) = 1
            Hel(4) = -2
            Hel(5) = 1
         end if


** The abbreviations - calculating vector products
         if (qk .eq. 1) then
            call hel_abbr_qq
         else if (qk .eq. 2) then
            call hel_abbr_qbarq
         else if (qk .eq. 3) then
            call hel_abbr_qqbar
         else if (qk .eq. 4) then
            call hel_abbr_qbarqbar
         end if
         
         
** Clearing the values of the tree and loop components
         Ctree(1) = (0D0,0D0)
         Cloop(1) = (0D0,0D0)

         CloopBuP(1) = (0D0,0D0)
         CloopBuN(1) = (0D0,0D0)
         CloopBlP(1) = (0D0,0D0)
         CloopBlN(1) = (0D0,0D0)
         CloopP1(1) = (0D0,0D0)
         CloopP2(1) = (0D0,0D0)
         
         CloopQct(1) = (0D0,0D0)
         CloopQlN(1) = (0D0,0D0)
         CloopQlP(1) = (0D0,0D0)
         CloopQuN(1) = (0D0,0D0)
         CloopQuP(1) = (0D0,0D0)


** Tree level
         call treeNC(mass2,gL,gR,S,T,U,T24,Ctree)

** Correcting the sign of the tree level component
         if (qk .eq. 2 .or. qk .eq. 3) then
            Ctree(1) = -Ctree(1)
         end if


** Calculating diagrams for the boxes attached to the lower quark line
         if (qk .eq. 1 .or. qk .eq. 2) then
            call boxesBlN_12(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBlN)
            call boxesBlP_12(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBlP)
         else if (qk .eq. 3 .or. qk .eq. 4) then
            call boxesBlN_34(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBlN)
            call boxesBlP_34(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBlP)
         end if

** Correcting the sign 
         if (NC .eq. 2 .or. NC .eq. 4) then
            CloopBlN(1) = -CloopBlN(1)
         end if
         if (qk .eq. 2 .or. qk .eq. 4) then
            CloopBlN(1) = -CloopBlN(1)
            CloopBlP(1) = -CloopBlP(1)
         end if
                             

** Calculating diagrams for the boxes attached to the upper quark line
         if (qk .eq. 1 .or. qk .eq. 3) then
            call boxesBuN_13(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBuN)
            call boxesBuP_13(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBuP)
         else if (qk .eq. 2 .or. qk .eq. 4) then
            call boxesBuN_24(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBuN)
            call boxesBuP_24(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopBuP)
         end if

** Correcting the sign 
         if (NC .eq. 3 .or. NC .eq. 4) then
            CloopBuN(1) = -CloopBuN(1)
         end if
         if (qk .eq. 3 .or. qk .eq. 4) then
            CloopBuN(1) = -CloopBuN(1)
            CloopBuP(1) = -CloopBuP(1)
         end if
         

** Calculating diagrams for the pentagons
         if (qk .eq. 1) then
            call pentP1_qq(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,CloopP1)
            if (cur .eq. 1 .or. cur .eq. 4) then
               call pentP2_qq14(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -              CloopP2)
            else if (cur .eq. 2 .or. cur .eq. 3) then
               call pentP2_qq23(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -              CloopP2)
            end if
         else if (qk .eq. 2) then
            call pentP1_qbarq(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopP1)
            if (cur .eq. 1 .or. cur .eq. 4) then
               call pentP2_qbarq14(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -              CloopP2)
            else if (cur .eq. 2 .or. cur .eq. 3) then
               call pentP2_qbarq23(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -              CloopP2)
            end if
         else if ( qk .eq. 3) then
            call pentP1_qqbar(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopP1)
            if (cur .eq. 1 .or. cur .eq. 4) then
               call pentP2_qqbar14(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -              CloopP2)
            else if (cur .eq. 2 .or. cur .eq. 3) then
               call pentP2_qqbar23(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -              CloopP2)
            end if
         else if (qk .eq. 4) then
            call pentP1_qbarqbar(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopP1)
            if (cur .eq. 1 .or. cur .eq. 4) then
               call pentP2_qbarqbar14(mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -              gR,CloopP2)
            else if (cur .eq. 2 .or. cur .eq. 3) then
               call pentP2_qbarqbar23(mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -              gR,CloopP2)
            end if
         end if
                          

** Diagram sets for the corrections to the qqV vertex (this is called iff 
** we're not using an effective qqV coupling)
         if (qqVmat) then
            call boxesQct(NC,mass2,gL,gR,q,dgL,dgR,S,T,T24,T14,S34,U,
     -           CloopQct)
            
            if (qk .eq. 1 .or. qk .eq. 2) then
               call boxesQlN_12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQlN)   
               call boxesQlP_12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQlP)
            else if (qk .eq. 3 .or. qk .eq. 4) then
               call boxesQlN_34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQlN)   
               call boxesQlP_34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQlP)
            end if
            
            if (qk .eq. 1 .or. qk .eq. 3) then     
               call boxesQuN_13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQuN)      
               call boxesQuP_13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQuP)
            else if (qk .eq. 2 .or. qk .eq. 4) then     
               call boxesQuN_24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQuN)      
               call boxesQuP_24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopQuP)
            end if
         end if
                                 
** Fixing the sign of the qqV loop components
         if (NC .eq. 2 .or. NC .eq. 4) then
            CloopQlN(1) = -CloopQlN(1)
         end if
         if (NC .eq. 3 .or. NC .eq. 4) then
            CloopQuN(1) = -CloopQuN(1)
         end if
                              
         if (qk .eq. 2 .or. qk .eq. 3) then
            CloopQct(1) = -CloopQct(1)
         end if
         
         if (qk .eq. 2 .or. qk .eq. 4) then
            CloopQlP(1) = -CloopQlP(1)
            CloopQlN(1) = -CloopQlN(1)
         end if
         
         if (qk .eq. 3 .or. qk .eq. 4) then
            CloopQuP(1) = -CloopQuP(1)
            CloopQuN(1) = -CloopQuN(1)
         end if
                          

** Setting individual loop components to zero for the purposes of debugging
*         CloopQct(1) = (0D0,0D0)
*         CloopQuN(1) = (0D0,0D0)
*         CloopQuP(1) = (0D0,0D0)
*         CloopQlN(1) = (0D0,0D0)
*         CloopQlP(1) = (0D0,0D0)
         
*         CloopBuP(1) = (0D0,0D0)
*         CloopBuN(1) = (0D0,0D0)
*         CloopBlP(1) = (0D0,0D0)
*         CloopBlN(1) = (0D0,0D0)
*         CloopP1(1) = (0D0,0D0)
*         CloopP2(1) = (0D0,0D0)
         

** The total loop component
         Cloop(1) = CloopBuP(1) + CloopBuN(1) + CloopBlP(1) +
     -        CloopBlN(1) + CloopP1(1) + CloopP2(1) + CloopQct(1) + 
     -        CloopQuP(1) + CloopQuN(1) + CloopQlP(1) + CloopQlN(1)


** Multiplying by SUSY factors, if necessary
         CtreeD1(1) = susyFAC(1)*Ctree(1)
         CtreeD2(1) = susyFAC(2)*Ctree(1)

         if (matFAC) then
*            SubMatEl = SubMatEl + 
*     -           sumup(MatSUN,susyFAC(1)*Ctree,susyFAC(2)*Ctree)  
            SubMatEl = SubMatEl + sumup(MatSUN,CtreeD1,CtreeD2)  
            Ctree(1) = susyFAC(1)*Ctree(1)
            Cloop(1) = susyFAC(2)*Cloop(1)
         else
            SubMatEl = SubMatEl + sumup(MatSUN,CtreeD1,CtreeD1)  
*            SubMatEl = SubMatEl + 
*     -           sumup(MatSUN,susyFAC(1)*Ctree,susyFAC(1)*Ctree)  
         end if
         

** The matrix elements squared
         SqMatEl(NC,1) = SqMatEl(NC,1) + sumup(MatSUN,Ctree,Ctree)
         
         SqMatEl(NC,2) = SqMatEl(NC,2) +2*sumup(MatSUN,Cloop,Ctree)
         

** Debugging messages 
         if (debug .gt. 1 .and. NC .eq. 1) then
            write(*,*)'Hel(i) =', Hel
            
            write(*,*)'Ctree =', Ctree(1)
            
            write(*,*)'CloopBuP =', CloopBuP(1) 
            write(*,*)'CloopBuN =', CloopBuP(1) + CloopBuN(1)
            write(*,*)'CloopBlP =', CloopBuP(1) + CloopBuN(1) + 
     -           CloopBlP(1)
            write(*,*)'CloopBlN =', CloopBuP(1) + CloopBuN(1) + 
     -           CloopBlP(1) + CloopBlN(1)
            write(*,*)'CloopP1 =', CloopBuP(1) + CloopBuN(1) + 
     -           CloopBlP(1) + CloopBlN(1) +  CloopP1(1)
            write(*,*)'CloopP2 =', CloopBuP(1) + CloopBuN(1) + 
     -           CloopBlP(1) + CloopBlN(1) +  CloopP1(1) + CloopP2(1)
            
            write(*,*)'CloopCT =', CloopQct(1) 
            write(*,*)'CloopQuP =', CloopQct(1) + CloopQuP(1)
            write(*,*)'CloopQuN =', CloopQct(1) + CloopQuP(1) +
     -           CloopQuN(1)
            write(*,*)'CloopQlP =', CloopQct(1) + CloopQuP(1) + 
     -           CloopQuN(1) + CloopQlP(1)
            write(*,*)'CloopQlN =', CloopQct(1) + CloopQuP(1) + 
     -          CloopQuN(1) + CloopQlP(1) + CloopQlN(1)
            
            write(*,*)'Cloop =', Cloop
            
            write(*,*)'Dsub =', Dsub
            write(*,*)'SqMatEl(NC,1) =', SqMatEl(NC,1)
            write(*,*)'SqMatEl(NC,2) =', SqMatEl(NC,2)
            
         end if
                          

** This will stop the code and output individual components if we get silly
** answers for the loop
         if (sqMatEl(NC,2) .eq. 1000D0) then
            write(*,*)'oops, NC, qk =', NC, qk    
            
            write(*,*)'Hel(i) =', Hel
            
            write(*,*)'Ctree =', Ctree(1)
            
            write(*,*)'CloopBuP =', CloopBuP(1) 
            write(*,*)'CloopBuN =', CloopBuP(1) + CloopBuN(1)
            write(*,*)'CloopBlP =', CloopBuP(1) +  CloopBuN(1) + 
     -           CloopBlP(1)
            write(*,*)'CloopBlN =', CloopBuP(1) + CloopBuN(1) +
     -           CloopBlP(1) + CloopBlN(1)
            write(*,*)'CloopP1 =', CloopBuP(1) + CloopBuN(1) + 
     -           CloopBlP(1) + CloopBlN(1) +  CloopP1(1)
            write(*,*)'CloopP2 =', CloopBuP(1) + CloopBuN(1) + 
     -           CloopBlP(1) + CloopBlN(1) +  CloopP1(1) + CloopP2(1)
            
            write(*,*)'CloopCT =', CloopQct(1) 
            write(*,*)'CloopQuP =', CloopQct(1) + CloopQuP(1)
            write(*,*)'CloopQuN =', CloopQct(1) + CloopQuP(1) +
     -           CloopQuN(1)
            write(*,*)'CloopQlP =', CloopQct(1) + CloopQuP(1) + 
     -           CloopQuN(1) + CloopQlP(1)
            write(*,*)'CloopQlN =', CloopQct(1) + CloopQuP(1) + 
     -           CloopQuN(1) + CloopQlP(1) + CloopQlN(1)
            
            write(*,*)'Cloop =', Cloop
            
            write(*,*)'Dsub =', Dsub
            write(*,*)'SqMatEl(NC,1) =', SqMatEl(NC,1)
            write(*,*)'SqMatEl(NC,2) =', SqMatEl(NC,2)
*            stop
            
            SqMatEl(NC,1) = 0D0
            SqMatEl(NC,2) = 0D0
         end if

      end do


** Debugging messages 
      if (debug .ne. 0) then
         write(*,*)'   '
         write(*,*)'NC, qk =', NC, qk
         write(*,*)'Dsub =', Dsub
         write(*,*)'SqMatEl(NC,1) =', SqMatEl(NC,1)
         write(*,*)'SqMatEl(NC,2) =', SqMatEl(NC,2)         
      end if
      

** Adding the subtraction term 
      SqMatEl(NC,2) = SqMatEl(NC,2) + (SubMatEl*Dsub)     


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the tree part of the calculation

      subroutine treeNC(mass2,gL,gR,S,T,U,T24,Ctree)

      implicit none

** Momentum invariants
      double precision S, T, U, T24

** Tree level component
      double complex Ctree(1)

** Mass^2 of external particles
      double precision mass2(5)

** Left / right coupling of particles
      double precision gL(4), gR(4)


** QED constant
      double precision Alfa, Alfa2

** Parameter and abbreviation declarations
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      Ctree(1) = Ctree(1) + 
     -     1/CW2*(-(8*Alfa*MW*Pi*(Abb4*EL*gL(1)*gL(3))) - 
     -     8*Alfa*MW*Pi*(Abb2*EL*gL(3)*gR(1)) - 
     -     8*Alfa*MW*Pi*(Abb3*EL*gL(1)*gR(3)) - 
     -     8*Alfa*MW*Pi*(Abb1*EL*gR(1)*gR(3)))/
     -     (SW*(-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - 
     -     S - T24 - U))


      end


******************************************************************************
******************************************************************************

*************************************************************************
** This set of subroutines calculates the spinor products for each     **
** quark pattern. They are stolen from the FormCalc output files       **
** abbr0_hel.F and abbr1_hel.F.                                        **
*************************************************************************
      
******************************************************************************
******************************************************************************

*** Spinor products for q + q -> q + H + q

      subroutine hel_abbr_qq

      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel

** Vector products
      double complex F18, F5, F9, F3, F4, F15, F6, F10, F8, F30, F20
      double complex F31, F34, F2, F23, F12, F25, F22, F19, F27, F48
      double complex F40, F11, F21, F47, F42, F7, F1, F16, F17, F32
      double complex F35, F29, F14, F13, F28, F39, F26, F41, F43
      double complex F51, F52, F33, F44, F45, F46, F49, F50, F24
      double complex F36, F37, F38


** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS, vBxS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS
      external vBxS


** Contains the declarations of the abbreviations
#include "abbr_NC.h"

      F5 = vSxS(DottedSpinor(3,1,7),Spinor(2,1,6))
      F3 = vSxS(DottedSpinor(3,1,6),Spinor(2,1,7))
      F4 = vSxS(DottedSpinor(5,1,7),Spinor(1,1,6))
      F6 = vSxS(DottedSpinor(5,1,6),Spinor(1,1,7))
      F8 = VSeS(DottedSpinor(3,1,7),DottedSpinor(5,1,7))
      F2 = VSeS(DottedSpinor(3,1,6),DottedSpinor(5,1,6))
      F7 = VSeS(Spinor(1,1,7),Spinor(2,1,7))
      F1 = VSeS(Spinor(1,1,6),Spinor(2,1,6))

      F30 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     Spinor(1,1,7)))
      F31 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(1,1,7)))
      F34 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(2,1,7)))
      F23 = vSxS(DottedSpinor(3,1,6),vBxS(k(2),
     -     Spinor(1,1,6)))
      F25 = vSxS(DottedSpinor(3,1,6),vBxS(k(4),
     -     Spinor(1,1,6)))
      F22 = vSxS(DottedSpinor(3,1,6),vBxS(k(4),
     -     Spinor(2,1,6)))
      F27 = vSxS(DottedSpinor(5,1,7),vVxS(k(1),
     -     Spinor(2,1,7)))
      F48 = vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(1,1,7)))
      F40 = vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(2,1,7)))
      F21 = vSxS(DottedSpinor(5,1,6),vBxS(k(1),
     -     Spinor(2,1,6)))
      F47 = vSxS(DottedSpinor(5,1,6),vBxS(k(3),
     -     Spinor(1,1,6)))
      F42 = vSxS(DottedSpinor(5,1,6),vBxS(k(3),
     -     Spinor(2,1,6)))
      F16 = vSxS(DottedSpinor(3,1,7),vVeS(k(1),
     -     DottedSpinor(5,1,6)))
      F17 = vSxS(DottedSpinor(3,1,7),vVeS(k(2),
     -     DottedSpinor(5,1,6)))
      F32 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     vBxS(k(4),Spinor(2,1,6))))
      F35 = vSxS(DottedSpinor(3,1,7),vVeS(k(3),
     -     DottedSpinor(5,1,6)))
      F29 = vSxS(DottedSpinor(3,1,7),vVeS(k(4),
     -     DottedSpinor(5,1,6)))
      F14 = vSxS(DottedSpinor(3,1,6),vBeS(k(1),
     -     DottedSpinor(5,1,7)))
      F13 = vSxS(DottedSpinor(3,1,6),vBeS(k(2),
     -     DottedSpinor(5,1,7)))
      F28 = vSxS(DottedSpinor(3,1,6),vBxS(k(2),
     -     vVxS(k(4),Spinor(2,1,7))))
      F39 = vSxS(DottedSpinor(3,1,6),vBeS(k(3),
     -     DottedSpinor(5,1,7)))
      F26 = vSxS(DottedSpinor(3,1,6),vBeS(k(4),
     -     DottedSpinor(5,1,7)))
      F41 = vSxS(DottedSpinor(5,1,7),vVxS(k(1),
     -     vBxS(k(3),Spinor(1,1,6))))
      F43 = vSxS(DottedSpinor(5,1,6),vBxS(k(1),
     -     vVxS(k(3),Spinor(1,1,7))))
      F51 = vSxS(DottedSpinor(3,1,7),
     -     vVxS(k(2),vBeS(k(1),DottedSpinor(5,1,7))))
      F52 = vSxS(DottedSpinor(3,1,7),
     -     vVxS(k(2),vBeS(k(3),DottedSpinor(5,1,7))))
      F33 = vSxS(DottedSpinor(3,1,7),
     -     vVxS(k(2),vBeS(k(4),DottedSpinor(5,1,7))))
      F44 = vSxS(DottedSpinor(3,1,7),
     -     vVxS(k(3),vBeS(k(1),DottedSpinor(5,1,7))))
      F45 = vSxS(DottedSpinor(3,1,7),
     -     vVxS(k(4),vBeS(k(1),DottedSpinor(5,1,7))))
      F46 = vSxS(DottedSpinor(3,1,7),
     -     vVxS(k(4),vBeS(k(3),DottedSpinor(5,1,7))))
      F49 = vSxS(DottedSpinor(3,1,6),
     -     vBxS(k(2),vVeS(k(1),DottedSpinor(5,1,6))))
      F50 = vSxS(DottedSpinor(3,1,6),
     -     vBxS(k(2),vVeS(k(3),DottedSpinor(5,1,6))))
      F24 = vSxS(DottedSpinor(3,1,6),
     -     vBxS(k(2),vVeS(k(4),DottedSpinor(5,1,6))))
      F36 = vSxS(DottedSpinor(3,1,6),
     -     vBxS(k(3),vVeS(k(1),DottedSpinor(5,1,6))))
      F37 = vSxS(DottedSpinor(3,1,6),
     -     vBxS(k(4),vVeS(k(1),DottedSpinor(5,1,6))))
      F38 = vSxS(DottedSpinor(3,1,6),
     -     vBxS(k(4),vVeS(k(3),DottedSpinor(5,1,6))))


      Abb1 = F1*F2
      Abb2 = F3*F4
      Abb3 = F5*F6
      Abb4 = F7*F8

      Abb44 = F27*F30
      Abb45 = F27*F31
      Abb46 = F33*F7
      Abb64 = F44*F7
      Abb65 = F45*F7
      Abb66 = F46*F7
      Abb75 = F30*F40
      Abb79 = F51*F7
      Abb80 = F52*F7
      Abb82 = F31*F40

      Abb40 = F21*F30
      Abb41 = F21*F31
      Abb42 = F32*F6
      Abb62 = F31*F42
      Abb63 = F43*F5
      Abb74 = F30*F42

      Abb32 = F23*F27
      Abb33 = F25*F27
      Abb34 = F28*F4
      Abb58 = F25*F40
      Abb59 = F3*F41
      Abb72 = F23*F40

      Abb22 = F21*F23
      Abb23 = F1*F24
      Abb24 = F21*F25
      Abb51 = F1*F36
      Abb52 = F1*F37
      Abb53 = F1*F38
      Abb69 = F23*F42
      Abb77 = F1*F49
      Abb78 = F1*F50
      Abb81 = F25*F42


** Abbreviation sums for the pentagons
      AbbSumP246 = Abb1 - 3*Abb2
      AbbSumP346 = Abb1 - 2*Abb2
      AbbSumP138 = Abb1 - Abb2
      AbbSumP343 = Abb1 - 1/2.D0*Abb2
      AbbSumP241 = Abb1 - 1/3.D0*Abb2
      AbbSumP9 = Abb1 + 1/4.D0*Abb2
      AbbSumP11 = Abb1 + Abb2
      AbbSumP20 = Abb1 + 4*Abb2
      AbbSumP342 = -2*Abb22 + Abb23
      AbbSumP46 = -Abb22 + Abb23
      AbbSumP252 = Abb1 - 2*Abb3
      AbbSumP240 = Abb1 - 1/2.D0*Abb3
      AbbSumP373 = Abb1 + 1/7.D0*Abb3
      AbbSumP445 = Abb1 + 1/6.D0*Abb3
      AbbSumP137 = Abb1 + 1/5.D0*Abb3
      AbbSumP277 = Abb1 + 1/4.D0*Abb3
      AbbSumP8 = Abb1 + 1/3.D0*Abb3
      AbbSumP325 = Abb1 + 1/2.D0*Abb3
      AbbSumP14 = Abb1 + Abb3
      AbbSumP330 = Abb1 + 2*Abb3
      AbbSumP29 = Abb1 + 3*Abb3
      AbbSumP282 = Abb1 + 4*Abb3
      AbbSumP150 = Abb1 + 5*Abb3
      AbbSumP450 = Abb1 + 6*Abb3
      AbbSumP379 = Abb1 + 7*Abb3
      AbbSumP71 = Abb32 - 1/2.D0*Abb33
      AbbSumP50 = Abb32 - Abb34
      AbbSumP347 = Abb32 - 1/2.D0*Abb34
      AbbSumP242 = -2*Abb2 + Abb4
      AbbSumP257 = Abb4 - 1/2.D0*Abb2
      AbbSumP381 = Abb4 + 1/7.D0*Abb2
      AbbSumP453 = Abb4 + 1/6.D0*Abb2
      AbbSumP155 = Abb4 + 1/5.D0*Abb2
      AbbSumP285 = Abb4 + 1/4.D0*Abb2
      AbbSumP37 = Abb4 + 1/3.D0*Abb2
      AbbSumP333 = Abb4 + 1/2.D0*Abb2
      AbbSumP23 = Abb2 + Abb4
      AbbSumP326 = 2*Abb2 + Abb4
      AbbSumP16 = 3*Abb2 + Abb4
      AbbSumP278 = 4*Abb2 + Abb4
      AbbSumP140 = 5*Abb2 + Abb4
      AbbSumP446 = 6*Abb2 + Abb4
      AbbSumP375 = 7*Abb2 + Abb4
      AbbSumP249 = -3*Abb3 + Abb4
      AbbSumP349 = -2*Abb3 + Abb4
      AbbSumP147 = -Abb3 + Abb4
      AbbSumP354 = Abb4 - 1/2.D0*Abb3
      AbbSumP256 = Abb4 - 1/3.D0*Abb3
      AbbSumP36 = Abb4 + 1/4.D0*Abb3
      AbbSumP32 = Abb3 + Abb4
      AbbSumP26 = 4*Abb3 + Abb4
      AbbSumP74 = Abb40 - 1/2.D0*Abb41
      AbbSumP53 = Abb40 - Abb42
      AbbSumP351 = Abb40 - 1/2.D0*Abb42
      AbbSumP353 = -2*Abb44 + Abb46
      AbbSumP55 = -Abb44 + Abb46
      AbbSumP189 = Abb22 + Abb51
      AbbSumP68 = -Abb22 + Abb52
      AbbSumP260 = Abb51 + Abb53
      AbbSumP428 = -Abb32 + Abb59
      AbbSumP169 = Abb32 + Abb59
      AbbSumP262 = Abb59 + 1/2.D0*Abb58
      AbbSumP430 = -Abb40 + Abb63
      AbbSumP171 = Abb40 + Abb63
      AbbSumP264 = Abb63 + 1/2.D0*Abb62
      AbbSumP193 = Abb44 + Abb64
      AbbSumP75 = -Abb44 + Abb65
      AbbSumP265 = Abb64 + Abb66
      AbbSumP84 = Abb51 + Abb69
      AbbSumP104 = Abb32 + Abb72
      AbbSumP88 = Abb59 + Abb72
      AbbSumP105 = Abb40 + Abb74
      AbbSumP92 = Abb63 + Abb74
      AbbSumP94 = Abb64 + Abb75
      AbbSumP111 = -Abb22 + Abb77
      AbbSumP210 = Abb78 - 1/2.D0*Abb69
      AbbSumP112 = -Abb44 + Abb79
      AbbSumP211 = Abb80 - 1/2.D0*Abb75
      AbbSumP434 = Abb53 - 1/2.D0*Abb81
      AbbSumP435 = Abb66 - 1/2.D0*Abb82
      AbbSumP7 = Abb24 + AbbSumP342
      AbbSumP21 = AbbSumP71 - 1/2.D0*Abb34
      AbbSumP30 = AbbSumP74 - 1/2.D0*Abb42
      AbbSumP35 = Abb45 + AbbSumP353
      AbbSumP276 = -Abb24 + AbbSumP189
      AbbSumP357 = AbbSumP68 - 1/2.D0*Abb24
      AbbSumP280 = -Abb33 + AbbSumP169
      AbbSumP283 = -Abb41 + AbbSumP171
      AbbSumP284 = -Abb45 + AbbSumP193
      AbbSumP362 = AbbSumP75 - 1/2.D0*Abb45
      AbbSumP239 = -Abb69 + AbbSumP46
      AbbSumP396 = Abb58 + Abb59 - Abb72
      AbbSumP301 = AbbSumP262 - 1/2.D0*Abb72
      AbbSumP448 = Abb59 + 1/3.D0*Abb58 - 1/3.D0*Abb72
      AbbSumP247 = Abb72 + AbbSumP50
      AbbSumP328 = -Abb58 + AbbSumP88
      AbbSumP399 = Abb62 + Abb63 - Abb74
      AbbSumP303 = AbbSumP264 - 1/2.D0*Abb74
      AbbSumP451 = Abb63 + 1/3.D0*Abb62 - 1/3.D0*Abb74
      AbbSumP253 = Abb74 + AbbSumP53
      AbbSumP331 = -Abb62 + AbbSumP92
      AbbSumP255 = -Abb75 + AbbSumP55
      AbbSumP103 = -Abb78 + AbbSumP111
      AbbSumP299 = Abb78 - AbbSumP260
      AbbSumP106 = -Abb80 + AbbSumP112
      AbbSumP304 = Abb80 - AbbSumP265
      AbbSumP324 = -Abb81 + AbbSumP84
      AbbSumP175 = AbbSumP260 - 1/2.D0*Abb81
      AbbSumP444 = Abb51 - 1/3.D0*Abb69 + 1/3.D0*Abb81
      AbbSumP394 = Abb51 - Abb69 + Abb81
      AbbSumP332 = -Abb82 + AbbSumP94
      AbbSumP178 = AbbSumP265 - 1/2.D0*Abb82
      AbbSumP452 = Abb64 - 1/3.D0*Abb75 + 1/3.D0*Abb82
      AbbSumP401 = Abb64 - Abb75 + Abb82
      AbbSumP59 = -Abb52 + Abb53 + AbbSumP189
      AbbSumP61 = AbbSumP262 + AbbSumP71
      AbbSumP64 = AbbSumP264 + AbbSumP74
      AbbSumP65 = -Abb65 + Abb66 + AbbSumP193
      AbbSumP409 = -Abb34 + AbbSumP396
      AbbSumP144 = -Abb58 + AbbSumP247
      AbbSumP413 = -Abb42 + AbbSumP399
      AbbSumP151 = -Abb62 + AbbSumP253
      AbbSumP384 = -Abb52 + Abb77 - 1/2.D0*Abb22 + 1/2.D0*Abb24
      AbbSumP204 = Abb77 - AbbSumP210 - 3/2.D0*Abb22
      AbbSumP385 = -Abb65 + Abb79 - 1/2.D0*Abb44 + 1/2.D0*Abb45
      AbbSumP205 = Abb79 - AbbSumP211 - 3/2.D0*Abb44
      AbbSumP460 = AbbSumP210 - AbbSumP434
      AbbSumP135 = Abb81 + AbbSumP239
      AbbSumP461 = AbbSumP211 - AbbSumP435
      AbbSumP153 = Abb82 + AbbSumP255
      AbbSumP292 = AbbSumP301 + 1/2.D0*Abb32 - 1/2.D0*Abb33
      AbbSumP293 = AbbSumP303 + 1/2.D0*Abb40 - 1/2.D0*Abb41
      AbbSumP254 = 0 + 1/2.D0*0 - 1/2.D0*0
      AbbSumP291 = -Abb52 + Abb77 - AbbSumP299
      AbbSumP294 = -Abb65 + Abb79 - AbbSumP304
      AbbSumP419 = Abb23 - Abb52 + AbbSumP434 + 1/2.D0*Abb24
      AbbSumP423 = Abb46 - Abb65 + AbbSumP435 + 1/2.D0*Abb45
      AbbSumP167 = -Abb52 + AbbSumP189 + AbbSumP434 + 1/2.D0*Abb24
      AbbSumP426 = Abb51 + AbbSumP357 - AbbSumP434
      AbbSumP172 = -Abb65 + AbbSumP193 + AbbSumP435 + 1/2.D0*Abb45
      AbbSumP431 = Abb64 + AbbSumP362 - AbbSumP435
      AbbSumP404 = Abb23 + Abb24 - Abb51 - 2*Abb52 + 2*Abb53 + 
     -   Abb69 - 2*Abb81
      AbbSumP416 = Abb45 + Abb46 - Abb64 - 2*Abb65 + 2*Abb66 + 
     -   Abb75 - 2*Abb82
      AbbSumP469 = AbbSumP384 - AbbSumP460
      AbbSumP471 = AbbSumP385 - AbbSumP461
      AbbSumP437 = Abb23 + Abb24 - 2*Abb52 + 2*Abb53 + Abb77 - 
     -   Abb81 - AbbSumP210 - 1/2.D0*Abb22
      AbbSumP463 = Abb53 + Abb69 - Abb78 - Abb81 + AbbSumP384 - 
     -     1/2.D0*Abb51
      AbbSumP455 = -Abb51 + AbbSumP469
      AbbSumP339 = AbbSumP175 - AbbSumP210 + AbbSumP384
      AbbSumP441 = Abb45 + Abb46 - 2*Abb65 + 2*Abb66 + Abb79 - 
     -   Abb82 - AbbSumP211 - 1/2.D0*Abb44
      AbbSumP467 = Abb66 + Abb75 - Abb80 - Abb82 + AbbSumP385 - 
     -     1/2.D0*Abb64
      AbbSumP456 = -Abb64 + AbbSumP471
      AbbSumP340 = AbbSumP178 - AbbSumP211 + AbbSumP385

** Abbreviation sums for boxes attached to the lower quark line
      AbbSumL30 = Abb40 - Abb74
      AbbSumL24 = Abb22 - Abb69
      AbbSumL70 = Abb1 - 4*Abb2
      AbbSumL76 = Abb1 - 1/2.D0*Abb2
      AbbSumL96 = Abb1 + 2*Abb2
      AbbSumL100 = Abb1 + 3*Abb2
      AbbSumL33 = Abb44 - Abb75
      AbbSumL28 = Abb32 - Abb72
      AbbSumL1 = Abb1 + Abb3
      AbbSumL46 = Abb41 - Abb62
      AbbSumL43 = Abb24 - Abb81
      AbbSumL48 = Abb45 - Abb82
      AbbSumL45 = Abb33 - Abb58
      AbbSumL2 = Abb2 + Abb4
      AbbSumL77 = -2*Abb3 + Abb4
      AbbSumL74 = Abb4 - 1/4.D0*Abb3
      AbbSumL102 = Abb4 + 1/3.D0*Abb3
      AbbSumL97 = Abb4 + 1/2.D0*Abb3
      AbbSumL4 = AbbSumL24 + AbbSumL30
      AbbSumL7 = AbbSumL28 + AbbSumL33
      AbbSumL55 = AbbSumL30 + AbbSumL46
      AbbSumL53 = AbbSumL24 + AbbSumL43
      AbbSumL11 = AbbSumL43 + AbbSumL46
      AbbSumL56 = AbbSumL33 + AbbSumL48
      AbbSumL87 = AbbSumL33 - 1/3.D0*Abb45 + 1/3.D0*Abb82
      AbbSumL114 = AbbSumL33 - 2/3.D0*Abb45 + 2/3.D0*Abb82
      AbbSumL75 = AbbSumL33 - AbbSumL48
      AbbSumL54 = AbbSumL28 + AbbSumL45
      AbbSumL13 = AbbSumL45 + AbbSumL48
      AbbSumL86 = AbbSumL28 - 1/3.D0*Abb33 + 1/3.D0*Abb58
      AbbSumL112 = AbbSumL28 - 2/3.D0*Abb33 + 2/3.D0*Abb58
      AbbSumL72 = AbbSumL28 - AbbSumL45
      AbbSumL66 = -2*Abb1 - 2*Abb3 + AbbSumL2
      AbbSumL64 = AbbSumL2 - 1/4.D0*Abb1 - 1/4.D0*Abb3
      AbbSumL93 = AbbSumL2 + 1/3.D0*Abb1 + 1/3.D0*Abb3
      AbbSumL90 = AbbSumL2 + 1/2.D0*Abb1 + 1/2.D0*Abb3
      AbbSumL16 = AbbSumL11 + AbbSumL4
      AbbSumL17 = AbbSumL13 + AbbSumL7
      AbbSumL80 = AbbSumL7 - 1/3.D0*Abb45 + 1/3.D0*Abb82 - 
     -     1/3.D0*Abb33 + 1/3.D0*Abb58
      AbbSumL107 = AbbSumL7 - 2/3.D0*Abb45 + 2/3.D0*Abb82 - 
     -     2/3.D0*Abb33 + 2/3.D0*Abb58
      AbbSumL65 = -AbbSumL13 + AbbSumL7


** Abbreviation sums for boxes attached to the upper line
      AbbSumU4 = Abb82 + Abb62
      AbbSumU3 = Abb1 + Abb2
      AbbSumU39 = -Abb82 + Abb45
      AbbSumU121 = Abb45 - 1/2.D0*Abb82
      AbbSumU43 = Abb45 + 1/2.D0*Abb82
      AbbSumU72 = 2*Abb82 + Abb45
      AbbSumU37 = -Abb62 + Abb41
      AbbSumU120 = Abb41 - 1/2.D0*Abb62
      AbbSumU41 = Abb41 + 1/2.D0*Abb62
      AbbSumU70 = 2*Abb62 + Abb41
      AbbSumU50 = Abb45 + Abb41
      AbbSumU77 = Abb33 + Abb24
      AbbSumU115 = Abb33 - Abb58
      AbbSumU119 = Abb33 - 1/2.D0*Abb58
      AbbSumU12 = Abb1 - Abb3
      AbbSumU126 = Abb1 + 1/4.D0*Abb3
      AbbSumU96 = Abb1 + 1/3.D0*Abb3
      AbbSumU60 = Abb1 + 2*Abb3
      AbbSumU29 = Abb1 + 7/3.D0*Abb3
      AbbSumU103 = Abb1 + 3*Abb3
      AbbSumU131 = Abb1 + 4*Abb3
      AbbSumU114 = Abb24 - Abb81
      AbbSumU117 = Abb24 - 1/2.D0*Abb81
      AbbSumU14 = -Abb2 + Abb4
      AbbSumU132 = Abb4 + 1/4.D0*Abb2
      AbbSumU104 = Abb4 + 1/3.D0*Abb2
      AbbSumU31 = Abb4 + 3/7.D0*Abb2
      AbbSumU62 = Abb4 + 1/2.D0*Abb2
      AbbSumU98 = 3*Abb2 + Abb4
      AbbSumU128 = 4*Abb2 + Abb4
      AbbSumU1 = Abb3 + Abb4
      AbbSumU23 = -AbbSumU4 + AbbSumU50
      AbbSumU87 = AbbSumU120 + AbbSumU121
      AbbSumU25 = AbbSumU41 + AbbSumU43
      AbbSumU58 = AbbSumU70 + AbbSumU72
      AbbSumU84 = -Abb58 - Abb81 + AbbSumU77
      AbbSumU86 = AbbSumU77 - 1/2.D0*Abb58 - 1/2.D0*Abb81
      AbbSumU5 = AbbSumU1 - AbbSumU3
      AbbSumU94 = Abb3 + AbbSumU132 + 1/4.D0*Abb1
      AbbSumU79 = Abb3 + AbbSumU104 + 1/3.D0*Abb1
      AbbSumU19 = Abb3 + AbbSumU31 + 3/7.D0*Abb1
      AbbSumU53 = Abb3 + AbbSumU62 + 1/2.D0*Abb1
      AbbSumU75 = 3*Abb1 + Abb3 + AbbSumU98
      AbbSumU92 = 4*Abb1 + Abb3 + AbbSumU128


** Abbreviation sums associated with corrections to the qqV vertex
      AbbSumQ35 = Abb1 - Abb2
      AbbSumQ60 = Abb1 - 1/3.D0*Abb2
      AbbSumQ9 = Abb1 + 1/3.D0*Abb3
      AbbSumQ1 = Abb1 + 1/2.D0*Abb3
      AbbSumQ25 = Abb1 + 5/3.D0*Abb3
      AbbSumQ5 = Abb1 + 2*Abb3
      AbbSumQ15 = Abb1 + 3*Abb3
      AbbSumQ17 = Abb4 + 1/3.D0*Abb2
      AbbSumQ7 = Abb4 + 1/2.D0*Abb2
      AbbSumQ26 = Abb4 + 3/5.D0*Abb2
      AbbSumQ3 = 2*Abb2 + Abb4
      AbbSumQ12 = 3*Abb2 + Abb4
      AbbSumQ61 = -3*Abb3 + Abb4
      AbbSumQ40 = -Abb3 + Abb4

      
      end


******************************************************************************
******************************************************************************

*** Spinor products for qbar + q -> qbar + H + q (qk = 2)

      subroutine hel_abbr_qbarq

      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel

** Vector products
      double complex F5, F11, F3, F17, F13, F4, F12, F6, F8, F15
      double complex F35, F30, F36, F2, F18, F25, F37, F27, F28, F20
      double complex F29, F39, F22, F10, F23, F24, F7, F1, F16, F9
      double complex F38, F34, F21, F14, F19, F33, F31, F32, F44
      double complex F45, F51, F52, F40, F46, F47, F48, F49, F50
      double complex F26, F41, F42, F43


** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS, vBxS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS
      external vBxS


** Contains the declarations of the abbreviations
#include "abbr_NC.h"
      F5 = Hel(1)*vSxS(DottedSpinor(1,-1,7),Spinor(2,1,6))
      F3 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),Spinor(2,1,7)))
      F4 = -(Hel(3)*vSxS(DottedSpinor(5,1,7),Spinor(3,-1,6)))
      F6 = Hel(3)*vSxS(DottedSpinor(5,1,6),Spinor(3,-1,7))
      F8 = Hel(1)*VSeS(DottedSpinor(1,-1,7),
     -     DottedSpinor(5,1,7))
      F2 = -(Hel(1)*VSeS(DottedSpinor(1,-1,6),
     -     DottedSpinor(5,1,6)))
      F7 = Hel(3)*VSeS(Spinor(3,-1,7),Spinor(2,1,7))
      F1 = -(Hel(3)*VSeS(Spinor(3,-1,6),Spinor(2,1,6)))

      F35 = Hel(1)*Hel(3)*
     -   vSxS(DottedSpinor(1,-1,7),vVxS(k(2),
     -     Spinor(3,-1,7)))
      F30 = Hel(1)*vSxS(DottedSpinor(1,-1,7),vVxS(k(4),
     -     Spinor(2,1,7)))
      F36 = Hel(1)*Hel(3)*
     -   vSxS(DottedSpinor(1,-1,7),vVxS(k(4),
     -     Spinor(3,-1,7)))
      F25 = Hel(1)*Hel(3)*
     -   vSxS(DottedSpinor(1,-1,6),vBxS(k(2),
     -     Spinor(3,-1,6)))
      F37 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(4),Spinor(2,1,6))))
      F27 = Hel(1)*Hel(3)*
     -   vSxS(DottedSpinor(1,-1,6),vBxS(k(4),
     -     Spinor(3,-1,6)))
      F28 = vSxS(DottedSpinor(5,1,7),vVxS(k(1),
     -     Spinor(2,1,7)))
      F20 = Hel(3)*vSxS(DottedSpinor(5,1,7),vVxS(k(1),
     -     Spinor(3,-1,7)))
      F29 = vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(2,1,7)))
      F39 = Hel(3)*vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(3,-1,7)))
      F22 = vSxS(DottedSpinor(5,1,6),vBxS(k(1),
     -     Spinor(2,1,6)))
      F10 = -(Hel(3)*vSxS(DottedSpinor(5,1,6),
     -      vBxS(k(1),Spinor(3,-1,6))))
      F23 = vSxS(DottedSpinor(5,1,6),vBxS(k(3),
     -     Spinor(2,1,6)))
      F24 = -(Hel(3)*vSxS(DottedSpinor(5,1,6),
     -      vBxS(k(3),Spinor(3,-1,6))))
      F9 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVeS(k(2),DottedSpinor(5,1,6)))
      F38 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVxS(k(2),vBxS(k(4),Spinor(2,1,6))))
      F34 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVeS(k(3),DottedSpinor(5,1,6)))
      F21 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVeS(k(4),DottedSpinor(5,1,6)))
      F19 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBeS(k(2),DottedSpinor(5,1,7))))
      F33 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(2),vVxS(k(4),Spinor(2,1,7)))))
      F31 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBeS(k(3),DottedSpinor(5,1,7))))
      F32 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBeS(k(4),DottedSpinor(5,1,7))))
      F44 = -(Hel(3)*vSxS(DottedSpinor(5,1,7),
     -      vVxS(k(1),vBxS(k(3),Spinor(3,-1,6)))))
      F45 = Hel(3)*vSxS(DottedSpinor(5,1,6),
     -    vBxS(k(1),vVxS(k(3),Spinor(3,-1,7))))
      F51 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVxS(k(2),vBeS(k(1),DottedSpinor(5,1,7))))
      F52 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVxS(k(2),vBeS(k(3),DottedSpinor(5,1,7))))
      F40 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVxS(k(2),vBeS(k(4),DottedSpinor(5,1,7))))
      F46 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVxS(k(3),vBeS(k(1),DottedSpinor(5,1,7))))
      F47 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVxS(k(4),vBeS(k(1),DottedSpinor(5,1,7))))
      F48 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -    vVxS(k(4),vBeS(k(3),DottedSpinor(5,1,7))))
      F49 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(2),vVeS(k(1),DottedSpinor(5,1,6)))))
      F50 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(2),vVeS(k(3),DottedSpinor(5,1,6)))))
      F26 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(2),vVeS(k(4),DottedSpinor(5,1,6)))))
      F41 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(3),vVeS(k(1),DottedSpinor(5,1,6)))))
      F42 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(4),vVeS(k(1),DottedSpinor(5,1,6)))))
      F43 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -      vBxS(k(4),vVeS(k(3),DottedSpinor(5,1,6)))))


      Abb1 = F1*F2
      Abb2 = F3*F4
      Abb3 = F5*F6
      Abb4 = F7*F8

      Abb32 = F25*F28
      Abb33 = F27*F28
      Abb51 = F33*F4
      Abb58 = F27*F29
      Abb59 = F3*F44
      Abb72 = F25*F29

      Abb8 = F1*F26
      Abb22 = F22*F25
      Abb24 = F22*F27
      Abb64 = F1*F41
      Abb65 = F1*F42
      Abb66 = F1*F43
      Abb69 = F23*F25
      Abb77 = F1*F49
      Abb78 = F1*F50
      Abb81 = F23*F27

      Abb44 = F28*F35
      Abb45 = F28*F36
      Abb47 = F47*F7
      Abb56 = F40*F7
      Abb55 = F46*F7
      Abb53 = F48*F7
      Abb75 = F29*F35
      Abb79 = F51*F7
      Abb80 = F52*F7
      Abb82 = F29*F36

      Abb18 = F45*F5
      Abb40 = F22*F35
      Abb41 = F22*F36
      Abb50 = F38*F6
      Abb62 = F23*F36
      Abb74 = F23*F35


      AbbSumP342 = Abb1 - 3*Abb2
      AbbSumP242 = Abb1 - 2*Abb2
      AbbSumP10 = Abb1 - Abb2
      AbbSumP237 = Abb1 - 1/2.D0*Abb2
      AbbSumP338 = Abb1 - 1/3.D0*Abb2
      AbbSumP139 = Abb1 + 1/4.D0*Abb2
      AbbSumP11 = Abb1 + Abb2
      AbbSumP146 = Abb1 + 4*Abb2
      AbbSumP45 = Abb69 + Abb8
      AbbSumP234 = 2*Abb69 + Abb8
      AbbSumP347 = Abb1 - 2*Abb3
      AbbSumP337 = Abb1 - 1/2.D0*Abb3
      AbbSumP287 = Abb1 + 1/7.D0*Abb3
      AbbSumP445 = Abb1 + 1/6.D0*Abb3
      AbbSumP236 = Abb1 + 1/4.D0*Abb3
      AbbSumP9 = Abb1 + 1/3.D0*Abb3
      AbbSumP317 = Abb1 + 1/2.D0*Abb3
      AbbSumP15 = Abb1 + Abb3
      AbbSumP322 = Abb1 + 2*Abb3
      AbbSumP29 = Abb1 + 3*Abb3
      AbbSumP248 = Abb1 + 4*Abb3
      AbbSumP450 = Abb1 + 6*Abb3
      AbbSumP293 = Abb1 + 7*Abb3
      AbbSumP205 = Abb32 + Abb72
      AbbSumP339 = -2*Abb2 + Abb4
      AbbSumP352 = Abb4 - 1/2.D0*Abb2
      AbbSumP295 = Abb4 + 1/7.D0*Abb2
      AbbSumP453 = Abb4 + 1/6.D0*Abb2
      AbbSumP485 = Abb4 + 1/5.D0*Abb2
      AbbSumP253 = Abb4 + 1/4.D0*Abb2
      AbbSumP37 = Abb4 + 1/3.D0*Abb2
      AbbSumP325 = Abb4 + 1/2.D0*Abb2
      AbbSumP23 = Abb2 + Abb4
      AbbSumP318 = 2*Abb2 + Abb4
      AbbSumP16 = 3*Abb2 + Abb4
      AbbSumP238 = 4*Abb2 + Abb4
      AbbSumP446 = 6*Abb2 + Abb4
      AbbSumP289 = 7*Abb2 + Abb4
      AbbSumP345 = -3*Abb3 + Abb4
      AbbSumP245 = -2*Abb3 + Abb4
      AbbSumP26 = -Abb3 + Abb4
      AbbSumP252 = Abb4 - 1/2.D0*Abb3
      AbbSumP351 = Abb4 - 1/3.D0*Abb3
      AbbSumP157 = Abb4 + 1/4.D0*Abb3
      AbbSumP32 = Abb3 + Abb4
      AbbSumP150 = 4*Abb3 + Abb4
      AbbSumP243 = Abb72 + 1/2.D0*Abb51
      AbbSumP50 = Abb72 + Abb51
      AbbSumP206 = Abb40 + Abb74
      AbbSumP249 = Abb74 + 1/2.D0*Abb50
      AbbSumP53 = Abb74 + Abb50
      AbbSumP54 = Abb75 + Abb56
      AbbSumP250 = 2*Abb75 + Abb56
      AbbSumP189 = Abb22 + Abb64
      AbbSumP82 = Abb69 + Abb64
      AbbSumP354 = Abb64 + Abb65
      AbbSumP177 = -Abb69 + Abb66
      AbbSumP191 = Abb32 + Abb59
      AbbSumP356 = Abb59 + 1/2.D0*Abb33
      AbbSumP61 = Abb72 + Abb59
      AbbSumP194 = Abb40 + Abb18
      AbbSumP64 = Abb74 + Abb18
      AbbSumP358 = Abb18 + 1/2.D0*Abb41
      AbbSumP196 = Abb44 + Abb55
      AbbSumP89 = Abb75 + Abb55
      AbbSumP359 = Abb55 + Abb47
      AbbSumP182 = -Abb75 + Abb53
      AbbSumP107 = Abb77 - 1/2.D0*Abb22
      AbbSumP212 = -Abb69 + Abb78
      AbbSumP108 = Abb79 - 1/2.D0*Abb44
      AbbSumP213 = -Abb75 + Abb80
      AbbSumP179 = Abb72 - 1/2.D0*Abb58
      AbbSumP181 = Abb74 - 1/2.D0*Abb62
      AbbSumP336 = Abb22 + AbbSumP45
      AbbSumP343 = Abb51 + AbbSumP205
      AbbSumP348 = Abb50 + AbbSumP206
      AbbSumP350 = Abb44 + AbbSumP54
      AbbSumP137 = -Abb81 + AbbSumP234
      AbbSumP274 = -Abb24 + AbbSumP189
      AbbSumP316 = -Abb81 + AbbSumP82
      AbbSumP444 = Abb64 - 1/3.D0*Abb69 + 1/3.D0*Abb81
      AbbSumP387 = -Abb69 + Abb81 + Abb64
      AbbSumP68 = AbbSumP354 - 1/2.D0*Abb24
      AbbSumP257 = AbbSumP177 - 1/2.D0*Abb81
      AbbSumP429 = Abb64 + AbbSumP177
      AbbSumP415 = Abb8 - Abb65 + Abb66
      AbbSumP276 = -Abb33 + AbbSumP191
      AbbSumP375 = AbbSumP356 - 1/2.D0*Abb32
      AbbSumP278 = -Abb41 + AbbSumP194
      AbbSumP377 = AbbSumP358 - 1/2.D0*Abb40
      AbbSumP155 = -Abb82 + AbbSumP250
      AbbSumP279 = -Abb45 + AbbSumP196
      AbbSumP324 = -Abb82 + AbbSumP89
      AbbSumP452 = Abb55 - 1/3.D0*Abb75 + 1/3.D0*Abb82
      AbbSumP397 = -Abb75 + Abb82 + Abb55
      AbbSumP71 = AbbSumP359 - 1/2.D0*Abb45
      AbbSumP262 = AbbSumP182 - 1/2.D0*Abb82
      AbbSumP434 = Abb55 + AbbSumP182
      AbbSumP420 = Abb56 - Abb47 + Abb53
      AbbSumP373 = Abb77 - AbbSumP354
      AbbSumP204 = Abb77 - AbbSumP212
      AbbSumP461 = -Abb64 - Abb66 + Abb78
      AbbSumP378 = Abb79 - AbbSumP359
      AbbSumP207 = Abb79 - AbbSumP213
      AbbSumP466 = -Abb55 - Abb53 + Abb80
      AbbSumP417 = Abb33 - 2*Abb51 - Abb58
      AbbSumP320 = -Abb58 + AbbSumP61
      AbbSumP147 = AbbSumP243 - 1/2.D0*Abb58
      AbbSumP448 = Abb59 - 1/3.D0*Abb72 + 1/3.D0*Abb58
      AbbSumP431 = Abb59 - AbbSumP179
      AbbSumP463 = Abb59 - 1/2.D0*Abb72 + 1/2.D0*Abb58
      AbbSumP390 = -Abb72 + Abb59 + Abb58
      AbbSumP419 = Abb41 - 2*Abb50 - Abb62
      AbbSumP323 = -Abb62 + AbbSumP64
      AbbSumP153 = AbbSumP249 - 1/2.D0*Abb62
      AbbSumP451 = Abb18 - 1/3.D0*Abb74 + 1/3.D0*Abb62
      AbbSumP433 = Abb18 - AbbSumP181
      AbbSumP465 = Abb18 - 1/2.D0*Abb74 + 1/2.D0*Abb62
      AbbSumP396 = -Abb74 + Abb18 + Abb62
      AbbSumP8 = -Abb24 + AbbSumP336
      AbbSumP20 = -Abb33 + AbbSumP343
      AbbSumP30 = -Abb41 + AbbSumP348
      AbbSumP36 = -Abb45 + AbbSumP350
      AbbSumP169 = Abb65 - Abb66 + AbbSumP82
      AbbSumP422 = -Abb65 + AbbSumP429
      AbbSumP174 = Abb47 - Abb53 + AbbSumP89
      AbbSumP427 = -Abb47 + AbbSumP434
      AbbSumP101 = -Abb78 + AbbSumP107 + 3/2.D0*Abb69
      AbbSumP477 = -Abb65 + Abb66 + Abb77 - Abb78
      AbbSumP298 = -Abb66 + Abb78 - 1/2.D0*Abb69 + 1/2.D0*Abb81
      AbbSumP102 = -Abb80 + AbbSumP108 + 3/2.D0*Abb75
      AbbSumP483 = -Abb47 + Abb53 + Abb79 - Abb80
      AbbSumP299 = -Abb53 + Abb80 - 1/2.D0*Abb75 + 1/2.D0*Abb82
      AbbSumP171 = AbbSumP179 + AbbSumP356
      AbbSumP424 = AbbSumP431 - 1/2.D0*Abb33
      AbbSumP479 = Abb32 - Abb33 - Abb72 + Abb58
      AbbSumP173 = AbbSumP64 + 1/2.D0*Abb41 - 1/2.D0*Abb62
      AbbSumP426 = AbbSumP433 - 1/2.D0*Abb41
      AbbSumP481 = Abb40 - Abb74 - Abb41 + Abb62
      AbbSumP331 = Abb66 - Abb78 + AbbSumP373
      AbbSumP455 = -Abb65 + Abb77 - AbbSumP461
      AbbSumP435 = Abb8 - 2*Abb65 + 2*Abb66 + Abb77 - Abb78
      AbbSumP334 = Abb53 - Abb80 + AbbSumP378
      AbbSumP458 = -Abb47 + Abb79 - AbbSumP466
      AbbSumP441 = Abb56 - 2*Abb47 + 2*Abb53 + Abb79 - Abb80
      AbbSumP332 = AbbSumP375 + 1/2.D0*Abb72 - 1/2.D0*Abb58
      AbbSumP456 = AbbSumP463 + 1/2.D0*Abb32 - 1/2.D0*Abb33
      AbbSumP438 = Abb32 - 2*Abb33 - Abb72 + 2*Abb51 + 2*Abb58
      AbbSumP470 = -2*Abb72 + 2*Abb58 + AbbSumP276
      AbbSumP333 = AbbSumP377 + 1/2.D0*Abb74 - 1/2.D0*Abb62
      AbbSumP457 = AbbSumP465 + 1/2.D0*Abb40 - 1/2.D0*Abb41
      AbbSumP440 = Abb40 - Abb74 - 2*Abb41 + 2*Abb50 + 2*Abb62
      AbbSumP473 = -2*Abb74 + 2*Abb62 + AbbSumP278
      AbbSumP58 = AbbSumP169 - 1/2.D0*Abb24 + 1/2.D0*Abb81
      AbbSumP65 = AbbSumP174 - 1/2.D0*Abb45 + 1/2.D0*Abb82
      AbbSumP405 = Abb51 + AbbSumP470
      AbbSumP410 = Abb50 + AbbSumP473
      AbbSumP272 = 0D0
      AbbSumP399 = -2*Abb69 + Abb8 + Abb81 - 2*Abb65 + 
     -     2*Abb66 + AbbSumP189
      AbbSumP412 = -2*Abb75 + Abb56 + Abb82 - 2*Abb47 + 
     -     2*Abb53 + AbbSumP196
      AbbSumP467 = AbbSumP477 - 1/2.D0*Abb69 + 1/2.D0*Abb81 + 1/2.D0*Abb64
      AbbSumP475 = AbbSumP483 - 1/2.D0*Abb75 + 1/2.D0*Abb82 + 1/2.D0*Abb55
      AbbSumP285 = AbbSumP107 - AbbSumP298 - AbbSumP68
      AbbSumP286 = AbbSumP108 - AbbSumP299 - AbbSumP71

** Abbreviation sums for boxes attached to the lower quark line
      AbbSumL30 = Abb40 - Abb74
      AbbSumL24 = Abb22 - Abb69
      AbbSumL70 = Abb1 - 4*Abb2
      AbbSumL76 = Abb1 - 1/2.D0*Abb2
      AbbSumL96 = Abb1 + 2*Abb2
      AbbSumL100 = Abb1 + 3*Abb2
      AbbSumL33 = Abb44 - Abb75
      AbbSumL28 = Abb32 - Abb72
      AbbSumL1 = Abb1 + Abb3
      AbbSumL46 = Abb41 - Abb62
      AbbSumL43 = Abb24 - Abb81
      AbbSumL48 = Abb45 - Abb82
      AbbSumL45 = Abb33 - Abb58
      AbbSumL2 = Abb2 + Abb4
      AbbSumL77 = -2*Abb3 + Abb4
      AbbSumL74 = Abb4 - 1/4.D0*Abb3
      AbbSumL102 = Abb4 + 1/3.D0*Abb3
      AbbSumL97 = Abb4 + 1/2.D0*Abb3
      AbbSumL4 = AbbSumL24 + AbbSumL30
      AbbSumL7 = AbbSumL28 + AbbSumL33
      AbbSumL55 = AbbSumL30 + AbbSumL46
      AbbSumL53 = AbbSumL24 + AbbSumL43
      AbbSumL11 = AbbSumL43 + AbbSumL46
      AbbSumL56 = AbbSumL33 + AbbSumL48
      AbbSumL87 = AbbSumL33 - 1/3.D0*Abb45 + 1/3.D0*Abb82
      AbbSumL114 = AbbSumL33 - 2/3.D0*Abb45 + 2/3.D0*Abb82
      AbbSumL75 = AbbSumL33 - AbbSumL48
      AbbSumL54 = AbbSumL28 + AbbSumL45
      AbbSumL13 = AbbSumL45 + AbbSumL48
      AbbSumL86 = AbbSumL28 - 1/3.D0*Abb33 + 1/3.D0*Abb58
      AbbSumL112 = AbbSumL28 - 2/3.D0*Abb33 + 2/3.D0*Abb58
      AbbSumL72 = AbbSumL28 - AbbSumL45
      AbbSumL66 = -2*Abb1 - 2*Abb3 + AbbSumL2
      AbbSumL64 = AbbSumL2 - 1/4.D0*Abb1 - 1/4.D0*Abb3
      AbbSumL93 = AbbSumL2 + 1/3.D0*Abb1 + 1/3.D0*Abb3
      AbbSumL90 = AbbSumL2 + 1/2.D0*Abb1 + 1/2.D0*Abb3
      AbbSumL16 = AbbSumL11 + AbbSumL4
      AbbSumL17 = AbbSumL13 + AbbSumL7
      AbbSumL80 = AbbSumL7 - 1/3.D0*Abb45 + 1/3.D0*Abb82 - 
     -     1/3.D0*Abb33 + 1/3.D0*Abb58
      AbbSumL107 = AbbSumL7 - 2/3.D0*Abb45 + 2/3.D0*Abb82 - 
     -     2/3.D0*Abb33 + 2/3.D0*Abb58
      AbbSumL65 = -AbbSumL13 + AbbSumL7

** Abbreviation sums for boxes attached to the upper quark line
      AbbSumU3 = Abb1 + Abb2
      AbbSumU37 = -Abb62 + Abb41
      AbbSumU123 = Abb41 - 1/2.D0*Abb62
      AbbSumU41 = Abb41 + 1/2.D0*Abb62
      AbbSumU72 = 2*Abb62 + Abb41
      AbbSumU50 = Abb45 + Abb41
      AbbSumU79 = Abb33 + Abb24
      AbbSumU118 = Abb33 - Abb58
      AbbSumU122 = Abb33 - 1/2.D0*Abb58
      AbbSumU11 = Abb1 - Abb3
      AbbSumU130 = Abb1 + 1/4.D0*Abb3
      AbbSumU99 = Abb1 + 1/3.D0*Abb3
      AbbSumU63 = Abb1 + 2*Abb3
      AbbSumU27 = Abb1 + 7/3.D0*Abb3
      AbbSumU105 = Abb1 + 3*Abb3
      AbbSumU133 = Abb1 + 4*Abb3
      AbbSumU117 = Abb24 - Abb81
      AbbSumU120 = Abb24 - 1/2.D0*Abb81
      AbbSumU12 = -Abb2 + Abb4
      AbbSumU135 = Abb4 + 1/4.D0*Abb2
      AbbSumU107 = Abb4 + 1/3.D0*Abb2
      AbbSumU29 = Abb4 + 3/7.D0*Abb2
      AbbSumU64 = Abb4 + 1/2.D0*Abb2
      AbbSumU100 = 3*Abb2 + Abb4
      AbbSumU131 = 4*Abb2 + Abb4
      AbbSumU1 = Abb3 + Abb4
      AbbSumU39 = Abb45 - Abb82
      AbbSumU124 = Abb45 - 1/2.D0*Abb82
      AbbSumU43 = Abb45 + 1/2.D0*Abb82
      AbbSumU4 = Abb62 + Abb82
      AbbSumU74 = Abb45 + 2*Abb82
      AbbSumU86 = -Abb58 - Abb81 + AbbSumU79
      AbbSumU88 = AbbSumU79 - 1/2.D0*Abb58 - 1/2.D0*Abb81
      AbbSumU5 = AbbSumU1 - AbbSumU3
      AbbSumU96 = Abb3 + AbbSumU135 + 1/4.D0*Abb1
      AbbSumU81 = Abb3 + AbbSumU107 + 1/3.D0*Abb1
      AbbSumU16 = Abb3 + AbbSumU29 + 3/7.D0*Abb1
      AbbSumU55 = Abb3 + AbbSumU64 + 1/2.D0*Abb1
      AbbSumU77 = 3*Abb1 + Abb3 + AbbSumU100
      AbbSumU94 = 4*Abb1 + Abb3 + AbbSumU131
      AbbSumU20 = AbbSumU37 + AbbSumU39
      AbbSumU89 = AbbSumU123 + AbbSumU124
      AbbSumU22 = AbbSumU41 + AbbSumU43
      AbbSumU60 = AbbSumU72 + AbbSumU74

** Abbreviation sums for the corrections to the qqV vertex
      AbbSumQ35 = Abb1 - Abb2
      AbbSumQ60 = Abb1 - 1/3.D0*Abb2
      AbbSumQ9 = Abb1 + 1/3.D0*Abb3
      AbbSumQ1 = Abb1 + 1/2.D0*Abb3
      AbbSumQ25 = Abb1 + 5/3.D0*Abb3
      AbbSumQ5 = Abb1 + 2*Abb3
      AbbSumQ15 = Abb1 + 3*Abb3
      AbbSumQ17 = Abb4 + 1/3.D0*Abb2
      AbbSumQ7 = Abb4 + 1/2.D0*Abb2
      AbbSumQ26 = Abb4 + 3/5.D0*Abb2
      AbbSumQ3 = 2*Abb2 + Abb4
      AbbSumQ12 = 3*Abb2 + Abb4
      AbbSumQ61 = -3*Abb3 + Abb4
      AbbSumQ40 = -Abb3 + Abb4


      end


******************************************************************************
******************************************************************************

*** Spinor products for q + qbar -> q + H + qbar (qk = 3)
      subroutine hel_abbr_qqbar
      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel

** Vector products
      double complex F3, F11, F5, F14, F16, F6, F12, F4, F8, F18
      double complex F27, F48, F41, F2, F15, F25, F47, F38, F29, F20
      double complex F32, F33, F23, F10, F22, F24, F7, F1, F9, F42
      double complex F13, F40, F21, F19, F43, F17, F39, F30, F31
      double complex F28, F51, F44, F45, F52, F46, F34, F49, F35
      double complex F36, F50, F37, F26


** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS, vBxS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS
      external vBxS


** Contains the declarations of the abbreviations
#include "abbr_NC.h"


      F3 = Hel(2)*vSxS(DottedSpinor(2,-1,7),Spinor(1,1,6))
      F5 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),Spinor(1,1,7)))
      F6 = -(Hel(5)*vSxS(DottedSpinor(3,1,7),Spinor(5,-1,6)))
      F4 = Hel(5)*vSxS(DottedSpinor(3,1,6),Spinor(5,-1,7))
      F8 = Hel(2)*VSeS(DottedSpinor(2,-1,7),
     -     DottedSpinor(3,1,7))
      F2 = -(Hel(2)*VSeS(DottedSpinor(2,-1,6),
     -     DottedSpinor(3,1,6)))
      F7 = Hel(5)*VSeS(Spinor(5,-1,7),Spinor(1,1,7))
      F1 = -(Hel(5)*VSeS(Spinor(5,-1,6),Spinor(1,1,6)))

      F27 = Hel(2)*Hel(5)*
     -   vSxS(DottedSpinor(2,-1,7),vVxS(k(1),
     -     Spinor(5,-1,7)))
      F48 = Hel(2)*vSxS(DottedSpinor(2,-1,7),vVxS(k(3),
     -     Spinor(1,1,7)))
      F41 = Hel(2)*Hel(5)*
     -   vSxS(DottedSpinor(2,-1,7),vVxS(k(3),
     -     Spinor(5,-1,7)))
      F25 = Hel(2)*Hel(5)*
     -   vSxS(DottedSpinor(2,-1,6),vBxS(k(1),
     -     Spinor(5,-1,6)))
      F47 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(3),Spinor(1,1,6))))
      F38 = Hel(2)*Hel(5)*
     -   vSxS(DottedSpinor(2,-1,6),vBxS(k(3),
     -     Spinor(5,-1,6)))
      F29 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     Spinor(1,1,7)))
      F20 = Hel(5)*vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     Spinor(5,-1,7)))
      F32 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(1,1,7)))
      F33 = Hel(5)*vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(5,-1,7)))
      F23 = vSxS(DottedSpinor(3,1,6),vBxS(k(2),
     -     Spinor(1,1,6)))
      F10 = -(Hel(5)*vSxS(DottedSpinor(3,1,6),
     -      vBxS(k(2),Spinor(5,-1,6))))
      F22 = vSxS(DottedSpinor(3,1,6),vBxS(k(4),
     -     Spinor(1,1,6)))
      F24 = -(Hel(5)*vSxS(DottedSpinor(3,1,6),
     -      vBxS(k(4),Spinor(5,-1,6))))
      F9 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVeS(k(1),DottedSpinor(3,1,6)))
      F42 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVxS(k(1),vBxS(k(3),Spinor(1,1,6))))
      F40 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVeS(k(3),DottedSpinor(3,1,6)))
      F21 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVeS(k(4),DottedSpinor(3,1,6)))
      F19 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBeS(k(1),DottedSpinor(3,1,7))))
      F43 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(1),vVxS(k(3),Spinor(1,1,7)))))
      F39 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBeS(k(3),DottedSpinor(3,1,7))))
      F30 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBeS(k(4),DottedSpinor(3,1,7))))
      F31 = -(Hel(5)*vSxS(DottedSpinor(3,1,7),
     -      vVxS(k(2),vBxS(k(4),Spinor(5,-1,6)))))
      F28 = Hel(5)*vSxS(DottedSpinor(3,1,6),
     -    vBxS(k(2),vVxS(k(4),Spinor(5,-1,7))))
      F51 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVxS(k(1),vBeS(k(2),DottedSpinor(3,1,7))))
      F44 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVxS(k(1),vBeS(k(3),DottedSpinor(3,1,7))))
      F45 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVxS(k(1),vBeS(k(4),DottedSpinor(3,1,7))))
      F52 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVxS(k(3),vBeS(k(2),DottedSpinor(3,1,7))))
      F46 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVxS(k(3),vBeS(k(4),DottedSpinor(3,1,7))))
      F34 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -    vVxS(k(4),vBeS(k(2),DottedSpinor(3,1,7))))
      F49 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(1),vVeS(k(2),DottedSpinor(3,1,6)))))
      F35 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(1),vVeS(k(3),DottedSpinor(3,1,6)))))
      F36 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(1),vVeS(k(4),DottedSpinor(3,1,6)))))
      F50 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(3),vVeS(k(2),DottedSpinor(3,1,6)))))
      F37 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(3),vVeS(k(4),DottedSpinor(3,1,6)))))
      F26 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -      vBxS(k(4),vVeS(k(2),DottedSpinor(3,1,6)))))


      Abb1 = F1*F2
      Abb2 = F3*F4
      Abb3 = F5*F6
      Abb4 = F7*F8

      Abb23 = F43*F6
      Abb38 = F25*F29
      Abb40 = F31*F5
      Abb41 = F25*F32
      Abb73 = F29*F38
      Abb62 = F32*F38

      Abb44 = F27*F29
      Abb45 = F27*F32
      Abb46 = F34*F7
      Abb63 = F44*F7
      Abb64 = F45*F7
      Abb65 = F46*F7
      Abb75 = F29*F41
      Abb79 = F51*F7
      Abb80 = F52*F7
      Abb82 = F32*F41

      Abb24 = F22*F25
      Abb25 = F23*F25
      Abb34 = F1*F26
      Abb49 = F1*F35
      Abb50 = F1*F36
      Abb51 = F1*F37
      Abb81 = F22*F38
      Abb70 = F23*F38
      Abb77 = F1*F49
      Abb78 = F1*F50

      Abb7 = F28*F3
      Abb33 = F22*F27
      Abb32 = F23*F27
      Abb58 = F22*F41
      Abb71 = F23*F41
      Abb55 = F4*F42


      AbbSumP160 = Abb1 - Abb2
      AbbSumP9 = Abb1 + 1/4.D0*Abb2
      AbbSumP265 = Abb1 + 1/3.D0*Abb2
      AbbSumP13 = Abb1 + Abb2
      AbbSumP271 = Abb1 + 3*Abb2
      AbbSumP20 = Abb1 + 4*Abb2
      AbbSumP355 = -2*Abb25 + Abb34
      AbbSumP50 = -Abb25 + Abb34
      AbbSumP342 = Abb1 - 2*Abb3
      AbbSumP96 = Abb1 - Abb3
      AbbSumP337 = Abb1 - 1/2.D0*Abb3
      AbbSumP312 = Abb1 + 1/5.D0*Abb3
      AbbSumP266 = Abb1 + 1/4.D0*Abb3
      AbbSumP10 = Abb1 + 1/3.D0*Abb3
      AbbSumP74 = Abb1 + 1/2.D0*Abb3
      AbbSumP395 = Abb1 + 2/3.D0*Abb3
      AbbSumP14 = Abb1 + Abb3
      AbbSumP401 = Abb1 + 3/2.D0*Abb3
      AbbSumP81 = Abb1 + 2*Abb3
      AbbSumP31 = Abb1 + 3*Abb3
      AbbSumP277 = Abb1 + 4*Abb3
      AbbSumP315 = Abb1 + 5*Abb3
      AbbSumP78 = -Abb33 + Abb32
      AbbSumP368 = Abb32 - 1/2.D0*Abb33
      AbbSumP358 = -2*Abb32 + Abb7
      AbbSumP54 = -Abb32 + Abb7
      AbbSumP82 = Abb38 - Abb41
      AbbSumP370 = Abb38 - 1/2.D0*Abb41
      AbbSumP339 = -2*Abb2 + Abb4
      AbbSumP98 = -Abb2 + Abb4
      AbbSumP346 = Abb4 - 1/2.D0*Abb2
      AbbSumP317 = Abb4 + 1/5.D0*Abb2
      AbbSumP281 = Abb4 + 1/4.D0*Abb2
      AbbSumP40 = Abb4 + 1/3.D0*Abb2
      AbbSumP85 = Abb4 + 1/2.D0*Abb2
      AbbSumP406 = Abb4 + 2/3.D0*Abb2
      AbbSumP25 = Abb2 + Abb4
      AbbSumP396 = Abb4 + 3/2.D0*Abb2
      AbbSumP76 = 2*Abb2 + Abb4
      AbbSumP17 = 3*Abb2 + Abb4
      AbbSumP268 = 4*Abb2 + Abb4
      AbbSumP314 = 5*Abb2 + Abb4
      AbbSumP168 = -Abb3 + Abb4
      AbbSumP41 = Abb4 + 1/4.D0*Abb3
      AbbSumP282 = Abb4 + 1/3.D0*Abb3
      AbbSumP33 = Abb3 + Abb4
      AbbSumP273 = 3*Abb3 + Abb4
      AbbSumP27 = 4*Abb3 + Abb4
      AbbSumP361 = -2*Abb38 + Abb40
      AbbSumP57 = -Abb38 + Abb40
      AbbSumP362 = -2*Abb44 + Abb46
      AbbSumP58 = -Abb44 + Abb46
      AbbSumP208 = -Abb25 + Abb49
      AbbSumP196 = Abb49 + Abb51
      AbbSumP212 = Abb32 - Abb55
      AbbSumP214 = Abb38 - Abb23
      AbbSumP215 = -Abb44 + Abb63
      AbbSumP202 = Abb63 + Abb65
      AbbSumP94 = Abb49 - Abb70
      AbbSumP227 = Abb32 - Abb71
      AbbSumP100 = -Abb55 + Abb71
      AbbSumP228 = Abb38 - Abb73
      AbbSumP105 = -Abb23 + Abb73
      AbbSumP108 = Abb63 - Abb75
      AbbSumP127 = Abb77 - 1/2.D0*Abb25
      AbbSumP226 = Abb77 - Abb78
      AbbSumP472 = -Abb51 + Abb78
      AbbSumP128 = Abb79 - 1/2.D0*Abb44
      AbbSumP229 = Abb79 - Abb80
      AbbSumP475 = -Abb65 + Abb80
      AbbSumP473 = Abb71 - Abb58
      AbbSumP199 = 2*Abb55 + Abb58
      AbbSumP474 = Abb73 - Abb62
      AbbSumP201 = 2*Abb23 + Abb62
      AbbSumP8 = Abb24 + AbbSumP355
      AbbSumP21 = Abb33 + AbbSumP358
      AbbSumP32 = Abb41 + AbbSumP361
      AbbSumP39 = Abb45 + AbbSumP362
      AbbSumP299 = Abb24 + AbbSumP208
      AbbSumP72 = Abb25 + Abb50 - 3/2.D0*Abb24
      AbbSumP365 = -Abb24 + Abb25 + Abb50
      AbbSumP285 = AbbSumP196 - 1/2.D0*Abb81
      AbbSumP65 = -Abb55 + AbbSumP78
      AbbSumP68 = -Abb23 + AbbSumP82
      AbbSumP303 = Abb45 + AbbSumP215
      AbbSumP83 = Abb44 + Abb64 - 3/2.D0*Abb45
      AbbSumP371 = Abb44 - Abb45 + Abb64
      AbbSumP289 = AbbSumP202 - 1/2.D0*Abb82
      AbbSumP264 = -Abb70 + AbbSumP50
      AbbSumP335 = Abb81 + AbbSumP94
      AbbSumP456 = Abb49 - 1/3.D0*Abb81 + 1/3.D0*Abb70
      AbbSumP393 = Abb49 - Abb81 + Abb70
      AbbSumP272 = -Abb71 + AbbSumP54
      AbbSumP278 = -Abb73 + AbbSumP57
      AbbSumP280 = -Abb75 + AbbSumP58
      AbbSumP344 = Abb82 + AbbSumP108
      AbbSumP462 = Abb63 - 1/3.D0*Abb82 + 1/3.D0*Abb75
      AbbSumP404 = -Abb82 + Abb63 + Abb75
      AbbSumP340 = -Abb58 + AbbSumP100
      AbbSumP398 = Abb55 + AbbSumP473
      AbbSumP459 = 3*Abb55 + AbbSumP473
      AbbSumP430 = Abb7 - 1/2.D0*Abb33 - 1/2.D0*Abb58
      AbbSumP343 = -Abb62 + AbbSumP105
      AbbSumP403 = Abb23 + AbbSumP474
      AbbSumP461 = 3*Abb23 + AbbSumP474
      AbbSumP432 = Abb40 - 1/2.D0*Abb41 - 1/2.D0*Abb62
      AbbSumP426 = -Abb24 + Abb34 + Abb50 - Abb51
      AbbSumP433 = -Abb45 + Abb46 + Abb64 - Abb65
      AbbSumP159 = Abb81 + AbbSumP264
      AbbSumP176 = Abb82 + AbbSumP280
      AbbSumP385 = Abb77 - AbbSumP365
      AbbSumP119 = -Abb78 + AbbSumP127 + 1/2.D0*Abb70
      AbbSumP386 = Abb79 - AbbSumP371
      AbbSumP120 = -Abb80 + AbbSumP128 + 1/2.D0*Abb75
      AbbSumP492 = AbbSumP473 + AbbSumP78
      AbbSumP191 = -Abb55 + AbbSumP368 - 1/2.D0*Abb58
      AbbSumP437 = Abb55 + AbbSumP368 - 1/2.D0*Abb58
      AbbSumP165 = -AbbSumP473 + AbbSumP54
      AbbSumP494 = AbbSumP474 + AbbSumP82
      AbbSumP193 = -Abb23 + AbbSumP370 - 1/2.D0*Abb62
      AbbSumP439 = Abb23 + AbbSumP370 - 1/2.D0*Abb62
      AbbSumP172 = -AbbSumP474 + AbbSumP57
      AbbSumP434 = Abb49 - Abb51 + AbbSumP365
      AbbSumP188 = -Abb50 + Abb51 + AbbSumP299
      AbbSumP440 = Abb63 - Abb65 + AbbSumP371
      AbbSumP194 = -Abb64 + Abb65 + AbbSumP303
      AbbSumP319 = -Abb78 + AbbSumP285 + 1/2.D0*Abb70
      AbbSumP320 = -Abb80 + AbbSumP289 + 1/2.D0*Abb75
      AbbSumP415 = -Abb33 + Abb7 + Abb55 + Abb71 - 2*Abb58
      AbbSumP481 = Abb55 + 2*Abb71 - 2*Abb58 + AbbSumP78
      AbbSumP451 = -Abb33 + Abb7 - Abb58 + 1/2.D0*Abb32 + 1/2.D0*Abb71
      AbbSumP350 = -2*Abb55 + AbbSumP492
      AbbSumP465 = 2*Abb55 + AbbSumP492
      AbbSumP421 = -Abb41 + Abb40 + Abb23 + Abb73 - 2*Abb62
      AbbSumP484 = Abb23 + 2*Abb73 - 2*Abb62 + AbbSumP82
      AbbSumP453 = -Abb41 + Abb40 - Abb62 + 1/2.D0*Abb38 + 1/2.D0*Abb73
      AbbSumP351 = -2*Abb23 + AbbSumP494
      AbbSumP466 = 2*Abb23 + AbbSumP494
      AbbSumP62 = -Abb50 + Abb51 + AbbSumP208 + 3/2.D0*Abb24 - 1/2.D0*Abb81
      AbbSumP69 = -Abb64 + Abb65 + AbbSumP215 + 3/2.D0*Abb45 - 1/2.D0*Abb82
      AbbSumP488 = Abb51 + AbbSumP226 - AbbSumP365
      AbbSumP495 = Abb65 + AbbSumP229 - AbbSumP371
      AbbSumP408 = -2*Abb24 + Abb34 + 2*Abb50 - 2*Abb51 + AbbSumP393
      AbbSumP423 = -2*Abb45 + Abb46 + 2*Abb64 - 2*Abb65 + AbbSumP404
      AbbSumP464 = AbbSumP226 - AbbSumP434
      AbbSumP349 = AbbSumP188 + AbbSumP226
      AbbSumP447 = 2*Abb24 - Abb25 - Abb34 - 2*Abb50 + 
     -     2*Abb51 + AbbSumP226
      AbbSumP467 = AbbSumP229 - AbbSumP440
      AbbSumP352 = AbbSumP194 + AbbSumP229
      AbbSumP454 = -Abb44 + 2*Abb45 - Abb46 - 2*Abb64 + 
     -     2*Abb65 + AbbSumP229
      AbbSumP476 = AbbSumP488 - 1/2.D0*Abb49 + 1/2.D0*Abb81 - 1/2.D0*Abb70
      AbbSumP309 = -Abb50 + AbbSumP226 + AbbSumP285 + 
     -     3/2.D0*Abb24 - 3/2.D0*Abb25 + 1/2.D0*Abb70
      AbbSumP486 = AbbSumP495 + 1/2.D0*Abb82 - 1/2.D0*Abb63 - 1/2.D0*Abb75
      AbbSumP310 = -Abb64 + AbbSumP229 + AbbSumP289 - 
     -     3/2.D0*Abb44 + 3/2.D0*Abb45 + 1/2.D0*Abb75

** Abbreviation sums for boxes attached to the lower quark line
      AbbSumL30 = Abb38 - Abb73
      AbbSumL23 = Abb25 - Abb70
      AbbSumL70 = Abb1 - 4*Abb2
      AbbSumL75 = Abb1 - 1/2.D0*Abb2
      AbbSumL95 = Abb1 + 2*Abb2
      AbbSumL100 = Abb1 + 3*Abb2
      AbbSumL33 = Abb44 - Abb75
      AbbSumL27 = Abb32 - Abb71
      AbbSumL1 = Abb1 + Abb3
      AbbSumL46 = Abb41 - Abb62
      AbbSumL43 = Abb24 - Abb81
      AbbSumL47 = Abb45 - Abb82
      AbbSumL44 = Abb33 - Abb58
      AbbSumL2 = Abb2 + Abb4
      AbbSumL76 = -2*Abb3 + Abb4
      AbbSumL72 = Abb4 - 1/4.D0*Abb3
      AbbSumL101 = Abb4 + 1/3.D0*Abb3
      AbbSumL97 = Abb4 + 1/2.D0*Abb3
      AbbSumL4 = AbbSumL23 + AbbSumL30
      AbbSumL7 = AbbSumL27 + AbbSumL33
      AbbSumL55 = AbbSumL30 + AbbSumL46
      AbbSumL53 = AbbSumL23 + AbbSumL43
      AbbSumL12 = AbbSumL43 + AbbSumL46
      AbbSumL56 = AbbSumL33 + AbbSumL47
      AbbSumL87 = AbbSumL33 - 1/3.D0*Abb45 + 1/3.D0*Abb82
      AbbSumL116 = AbbSumL33 - 2/3.D0*Abb45 + 2/3.D0*Abb82
      AbbSumL74 = AbbSumL33 - AbbSumL47
      AbbSumL54 = AbbSumL27 + AbbSumL44
      AbbSumL13 = AbbSumL44 + AbbSumL47
      AbbSumL86 = AbbSumL27 - 1/3.D0*Abb33 + 1/3.D0*Abb58
      AbbSumL114 = AbbSumL27 - 2/3.D0*Abb33 + 2/3.D0*Abb58
      AbbSumL71 = AbbSumL27 - AbbSumL44
      AbbSumL66 = -2*Abb1 - 2*Abb3 + AbbSumL2
      AbbSumL63 = AbbSumL2 - 1/4.D0*Abb1 - 1/4.D0*Abb3
      AbbSumL92 = AbbSumL2 + 1/3.D0*Abb1 + 1/3.D0*Abb3
      AbbSumL90 = AbbSumL2 + 1/2.D0*Abb1 + 1/2.D0*Abb3
      AbbSumL16 = AbbSumL12 + AbbSumL4
      AbbSumL17 = AbbSumL13 + AbbSumL7
      AbbSumL79 = AbbSumL7 - 1/3.D0*Abb45 + 1/3.D0*Abb82 - 
     -     1/3.D0*Abb33 + 1/3.D0*Abb58
      AbbSumL108 = AbbSumL7 - 2/3.D0*Abb45 + 2/3.D0*Abb82 - 
     -     2/3.D0*Abb33 + 2/3.D0*Abb58
      AbbSumL65 = -AbbSumL13 + AbbSumL7

** Abbreviation sums for boxes attached to the upper quark line
      AbbSumU4 = Abb82 + Abb62
      AbbSumU3 = Abb1 + Abb2
      AbbSumU39 = -Abb82 + Abb45
      AbbSumU121 = Abb45 - 1/2.D0*Abb82
      AbbSumU43 = Abb45 + 1/2.D0*Abb82
      AbbSumU72 = 2*Abb82 + Abb45
      AbbSumU37 = -Abb62 + Abb41
      AbbSumU120 = Abb41 - 1/2.D0*Abb62
      AbbSumU41 = Abb41 + 1/2.D0*Abb62
      AbbSumU70 = 2*Abb62 + Abb41
      AbbSumU50 = Abb45 + Abb41
      AbbSumU77 = Abb33 + Abb24
      AbbSumU115 = Abb33 - Abb58
      AbbSumU119 = Abb33 - 1/2.D0*Abb58
      AbbSumU12 = Abb1 - Abb3
      AbbSumU126 = Abb1 + 1/4.D0*Abb3
      AbbSumU96 = Abb1 + 1/3.D0*Abb3
      AbbSumU60 = Abb1 + 2*Abb3
      AbbSumU29 = Abb1 + 7/3.D0*Abb3
      AbbSumU103 = Abb1 + 3*Abb3
      AbbSumU131 = Abb1 + 4*Abb3
      AbbSumU114 = Abb24 - Abb81
      AbbSumU117 = Abb24 - 1/2.D0*Abb81
      AbbSumU14 = -Abb2 + Abb4
      AbbSumU132 = Abb4 + 1/4.D0*Abb2
      AbbSumU104 = Abb4 + 1/3.D0*Abb2
      AbbSumU31 = Abb4 + 3/7.D0*Abb2
      AbbSumU62 = Abb4 + 1/2.D0*Abb2
      AbbSumU98 = 3*Abb2 + Abb4
      AbbSumU128 = 4*Abb2 + Abb4
      AbbSumU1 = Abb3 + Abb4
      AbbSumU23 = -AbbSumU4 + AbbSumU50
      AbbSumU87 = AbbSumU120 + AbbSumU121
      AbbSumU25 = AbbSumU41 + AbbSumU43
      AbbSumU58 = AbbSumU70 + AbbSumU72
      AbbSumU84 = -Abb58 - Abb81 + AbbSumU77
      AbbSumU86 = AbbSumU77 - 1/2.D0*Abb58 - 1/2.D0*Abb81
      AbbSumU5 = AbbSumU1 - AbbSumU3
      AbbSumU94 = Abb3 + AbbSumU132 + 1/4.D0*Abb1
      AbbSumU79 = Abb3 + AbbSumU104 + 1/3.D0*Abb1
      AbbSumU19 = Abb3 + AbbSumU31 + 3/7.D0*Abb1
      AbbSumU53 = Abb3 + AbbSumU62 + 1/2.D0*Abb1
      AbbSumU75 = 3*Abb1 + Abb3 + AbbSumU98
      AbbSumU92 = 4*Abb1 + Abb3 + AbbSumU128

** Abbreviation sums for the corrections to the qqV vertex
      AbbSumQ35 = Abb1 - Abb2
      AbbSumQ60 = Abb1 - 1/3.D0*Abb2
      AbbSumQ9 = Abb1 + 1/3.D0*Abb3
      AbbSumQ1 = Abb1 + 1/2.D0*Abb3
      AbbSumQ25 = Abb1 + 5/3.D0*Abb3
      AbbSumQ5 = Abb1 + 2*Abb3
      AbbSumQ15 = Abb1 + 3*Abb3
      AbbSumQ17 = Abb4 + 1/3.D0*Abb2
      AbbSumQ7 = Abb4 + 1/2.D0*Abb2
      AbbSumQ26 = Abb4 + 3/5.D0*Abb2
      AbbSumQ3 = 2*Abb2 + Abb4
      AbbSumQ12 = 3*Abb2 + Abb4
      AbbSumQ61 = -3*Abb3 + Abb4
      AbbSumQ40 = -Abb3 + Abb4

      end


******************************************************************************
******************************************************************************

*** Spinor products for qbar + qbar -> qbar + H + qbar (qk = 4)

      subroutine hel_abbr_qbarqbar

      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel

** Vector products
      double complex F11, F5, F16, F3, F4, F17, F6, F12, F8, F35
      double complex F14, F37, F31, F2, F26, F18, F27, F36, F15, F29
      double complex F34, F30, F13, F22, F32, F23, F7, F1, F19, F9
      double complex F38, F39, F21, F10, F20, F33, F24, F25, F44
      double complex F45, F51, F52, F40, F46, F47, F48, F49, F50
      double complex F28, F41, F42, F43


** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS, vBxS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS
      external vBxS


** Contains the declarations of the abbreviations
#include "abbr_NC.h"


      F5 = -(Hel(1)*Hel(5)*vSxS(DottedSpinor(1,-1,7),
     -     Spinor(5,-1,6)))
      F3 = -(Hel(1)*Hel(5)*vSxS(DottedSpinor(1,-1,6),
     -     Spinor(5,-1,7)))
      F4 = -(Hel(2)*Hel(3)*vSxS(DottedSpinor(2,-1,7),
     -     Spinor(3,-1,6)))
      F6 = -(Hel(2)*Hel(3)*vSxS(DottedSpinor(2,-1,6),
     -     Spinor(3,-1,7)))
      F8 = Hel(1)*Hel(2)*
     -     VSeS(DottedSpinor(1,-1,7),DottedSpinor(2,-1,7))
      F2 = Hel(1)*Hel(2)*
     -     VSeS(DottedSpinor(1,-1,6),DottedSpinor(2,-1,6))
      F7 = Hel(3)*Hel(5)*VSeS(Spinor(3,-1,7),Spinor(5,-1,7))
      F1 = Hel(3)*Hel(5)*VSeS(Spinor(3,-1,6),Spinor(5,-1,6))

      F35 = Hel(1)*Hel(3)*
     -     vSxS(DottedSpinor(1,-1,7),vVxS(k(2),Spinor(3,-1,7)))
      F14 = Hel(1)*Hel(5)*
     -     vSxS(DottedSpinor(1,-1,7),vVxS(k(2),Spinor(5,-1,7)))
      F37 = Hel(1)*Hel(3)*
     -     vSxS(DottedSpinor(1,-1,7),vVxS(k(4),Spinor(3,-1,7)))
      F31 = Hel(1)*Hel(5)*
     -     vSxS(DottedSpinor(1,-1,7),vVxS(k(4),Spinor(5,-1,7)))
      F26 = Hel(1)*Hel(3)*
     -     vSxS(DottedSpinor(1,-1,6),vBxS(k(2),Spinor(3,-1,6)))
      F18 = Hel(1)*Hel(5)*
     -     vSxS(DottedSpinor(1,-1,6),vBxS(k(2),
     -     Spinor(5,-1,6)))
      F27 = Hel(1)*Hel(3)*
     -     vSxS(DottedSpinor(1,-1,6),vBxS(k(4),
     -     Spinor(3,-1,6)))
      F36 = Hel(1)*Hel(5)*
     -     vSxS(DottedSpinor(1,-1,6),vBxS(k(4),
     -     Spinor(5,-1,6)))
      F15 = Hel(2)*Hel(3)*
     -     vSxS(DottedSpinor(2,-1,7),vVxS(k(1),
     -     Spinor(3,-1,7)))
      F29 = Hel(2)*Hel(5)*
     -     vSxS(DottedSpinor(2,-1,7),vVxS(k(1),
     -     Spinor(5,-1,7)))
      F34 = Hel(2)*Hel(3)*
     -     vSxS(DottedSpinor(2,-1,7),vVxS(k(3),
     -     Spinor(3,-1,7)))
      F30 = Hel(2)*Hel(5)*
     -     vSxS(DottedSpinor(2,-1,7),vVxS(k(3),
     -     Spinor(5,-1,7)))
      F13 = Hel(2)*Hel(3)*
     -     vSxS(DottedSpinor(2,-1,6),vBxS(k(1),
     -     Spinor(3,-1,6)))
      F22 = Hel(2)*Hel(5)*
     -     vSxS(DottedSpinor(2,-1,6),vBxS(k(1),
     -     Spinor(5,-1,6)))
      F32 = Hel(2)*Hel(3)*
     -     vSxS(DottedSpinor(2,-1,6),vBxS(k(3),
     -     Spinor(3,-1,6)))
      F23 = Hel(2)*Hel(5)*
     -     vSxS(DottedSpinor(2,-1,6),vBxS(k(3),
     -     Spinor(5,-1,6)))
      F38 = -(Hel(1)*Hel(5)*
     -     vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBxS(k(4),Spinor(5,-1,6)))))
      F39 = -(Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),vVeS(k(3),
     -     DottedSpinor(2,-1,6))))
      F21 = -(Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),vVeS(k(4),
     -     DottedSpinor(2,-1,6))))
      F33 = -(Hel(1)*Hel(5)*
     -     vSxS(DottedSpinor(1,-1,6),
     -     vBxS(k(2),vVxS(k(4),Spinor(5,-1,7)))))
      F24 = -(Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),vBeS(k(3),
     -     DottedSpinor(2,-1,7))))
      F25 = -(Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),vBeS(k(4),
     -     DottedSpinor(2,-1,7))))
      F44 = -(Hel(2)*Hel(3)*
     -     vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(1),vBxS(k(3),Spinor(3,-1,6)))))
      F45 = -(Hel(2)*Hel(3)*
     -     vSxS(DottedSpinor(2,-1,6),
     -     vBxS(k(1),vVxS(k(3),Spinor(3,-1,7)))))
      F51 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(1),DottedSpinor(2,-1,7))))
      F52 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(3),DottedSpinor(2,-1,7))))
      F40 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(4),DottedSpinor(2,-1,7))))
      F46 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(3),vBeS(k(1),DottedSpinor(2,-1,7))))
      F47 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),vBeS(k(1),DottedSpinor(2,-1,7))))
      F48 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),vBeS(k(3),DottedSpinor(2,-1,7))))
      F49 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),
     -     vBxS(k(2),vVeS(k(1),DottedSpinor(2,-1,6))))
      F50 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),
     -     vBxS(k(2),vVeS(k(3),DottedSpinor(2,-1,6))))
      F28 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),
     -     vBxS(k(2),vVeS(k(4),DottedSpinor(2,-1,6))))
      F41 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),
     -     vBxS(k(3),vVeS(k(1),DottedSpinor(2,-1,6))))
      F42 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),
     -     vBxS(k(4),vVeS(k(1),DottedSpinor(2,-1,6))))
      F43 = Hel(1)*Hel(2)*
     -     vSxS(DottedSpinor(1,-1,6),
     -     vBxS(k(4),vVeS(k(3),DottedSpinor(2,-1,6))))


      Abb1 = F1*F2
      Abb2 = F3*F4
      Abb3 = F5*F6
      Abb4 = F7*F8

      Abb18 = F23*F26
      Abb24 = F22*F27
      Abb29 = F1*F28
      Abb59 = F22*F26
      Abb63 = F23*F27
      Abb64 = F1*F41
      Abb65 = F1*F42
      Abb66 = F1*F43
      Abb79 = F1*F49
      Abb80 = F1*F50

      Abb33 = F27*F29
      Abb39 = F26*F30
      Abb40 = F33*F4
      Abb50 = F26*F29
      Abb69 = F27*F30
      Abb70 = F3*F44

      Abb20 = F38*F6
      Abb41 = F22*F37
      Abb46 = F23*F35
      Abb55 = F22*F35
      Abb62 = F23*F37
      Abb72 = F45*F5

      Abb32 = F30*F35
      Abb45 = F29*F37
      Abb53 = F29*F35
      Abb56 = F40*F7
      Abb73 = F52*F7
      Abb74 = F46*F7
      Abb75 = F47*F7
      Abb76 = F48*F7
      Abb81 = F51*F7
      Abb82 = F30*F37


** Abbreviations for the pentagons
      AbbSumP421 = Abb1 - 2*Abb2
      AbbSumP6 = Abb1 - Abb2
      AbbSumP415 = Abb1 - 1/2.D0*Abb2
      AbbSumP153 = Abb1 + 1/4.D0*Abb2
      AbbSumP351 = Abb1 + 1/3.D0*Abb2
      AbbSumP9 = Abb1 + Abb2
      AbbSumP355 = Abb1 + 3*Abb2
      AbbSumP160 = Abb1 + 4*Abb2
      AbbSumP43 = Abb18 + Abb29
      AbbSumP271 = 2*Abb18 + Abb29
      AbbSumP72 = Abb1 - 3*Abb3
      AbbSumP337 = Abb1 - 2*Abb3
      AbbSumP86 = Abb1 - Abb3
      AbbSumP333 = Abb1 - 1/2.D0*Abb3
      AbbSumP67 = Abb1 - 1/3.D0*Abb3
      AbbSumP152 = Abb1 + 1/5.D0*Abb3
      AbbSumP350 = Abb1 + 1/4.D0*Abb3
      AbbSumP478 = Abb1 + 2/7.D0*Abb3
      AbbSumP5 = Abb1 + 1/3.D0*Abb3
      AbbSumP201 = Abb1 + 1/2.D0*Abb3
      AbbSumP396 = Abb1 + 2/3.D0*Abb3
      AbbSumP11 = Abb1 + Abb3
      AbbSumP406 = Abb1 + 3/2.D0*Abb3
      AbbSumP208 = Abb1 + 2*Abb3
      AbbSumP26 = Abb1 + 3*Abb3
      AbbSumP483 = Abb1 + 7/2.D0*Abb3
      AbbSumP360 = Abb1 + 4*Abb3
      AbbSumP169 = Abb1 + 5*Abb3
      AbbSumP108 = Abb50 - Abb39
      AbbSumP68 = -3*Abb2 + Abb4
      AbbSumP334 = -2*Abb2 + Abb4
      AbbSumP88 = -Abb2 + Abb4
      AbbSumP339 = Abb4 - 1/2.D0*Abb2
      AbbSumP75 = Abb4 - 1/3.D0*Abb2
      AbbSumP177 = Abb4 + 1/5.D0*Abb2
      AbbSumP365 = Abb4 + 1/4.D0*Abb2
      AbbSumP484 = Abb4 + 2/7.D0*Abb2
      AbbSumP35 = Abb4 + 1/3.D0*Abb2
      AbbSumP512 = Abb4 + 2/5.D0*Abb2
      AbbSumP211 = Abb4 + 1/2.D0*Abb2
      AbbSumP409 = Abb4 + 2/3.D0*Abb2
      AbbSumP19 = Abb2 + Abb4
      AbbSumP398 = Abb4 + 3/2.D0*Abb2
      AbbSumP202 = 2*Abb2 + Abb4
      AbbSumP12 = 3*Abb2 + Abb4
      AbbSumP480 = Abb4 + 7/2.D0*Abb2
      AbbSumP352 = 4*Abb2 + Abb4
      AbbSumP156 = 5*Abb2 + Abb4
      AbbSumP426 = -2*Abb3 + Abb4
      AbbSumP23 = -Abb3 + Abb4
      AbbSumP433 = Abb4 - 1/2.D0*Abb3
      AbbSumP176 = Abb4 + 1/4.D0*Abb3
      AbbSumP364 = Abb4 + 1/3.D0*Abb3
      AbbSumP28 = Abb3 + Abb4
      AbbSumP358 = 3*Abb3 + Abb4
      AbbSumP166 = 4*Abb3 + Abb4
      AbbSumP276 = Abb39 + 1/2.D0*Abb40
      AbbSumP48 = Abb39 + Abb40
      AbbSumP109 = Abb55 - Abb46
      AbbSumP280 = Abb46 + 1/2.D0*Abb20
      AbbSumP51 = Abb46 + Abb20
      AbbSumP52 = Abb32 + Abb56
      AbbSumP281 = 2*Abb32 + Abb56
      AbbSumP219 = -Abb59 + Abb64
      AbbSumP85 = -Abb18 + Abb64
      AbbSumP66 = Abb64 + Abb65
      AbbSumP205 = Abb39 - Abb69
      AbbSumP288 = Abb39 - 1/2.D0*Abb69
      AbbSumP441 = Abb40 + Abb69
      AbbSumP221 = -Abb50 + Abb70
      AbbSumP70 = Abb70 + 1/2.D0*Abb33
      AbbSumP90 = -Abb39 + Abb70
      AbbSumP209 = Abb46 - Abb62
      AbbSumP290 = Abb46 - 1/2.D0*Abb62
      AbbSumP443 = Abb20 + Abb62
      AbbSumP224 = -Abb55 + Abb72
      AbbSumP94 = -Abb46 + Abb72
      AbbSumP73 = Abb72 + 1/2.D0*Abb41
      AbbSumP226 = -Abb53 + Abb74
      AbbSumP96 = -Abb32 + Abb74
      AbbSumP74 = Abb74 + Abb75
      AbbSumP107 = Abb79 - Abb80
      AbbSumP244 = Abb80 - 1/2.D0*Abb18
      AbbSumP110 = Abb81 - Abb73
      AbbSumP245 = Abb73 - 1/2.D0*Abb32
      AbbSumP349 = Abb59 + AbbSumP43
      AbbSumP356 = Abb50 + AbbSumP48
      AbbSumP361 = Abb55 + AbbSumP51
      AbbSumP363 = Abb53 + AbbSumP52
      AbbSumP150 = -Abb63 + AbbSumP271
      AbbSumP302 = Abb24 + AbbSumP219
      AbbSumP395 = Abb18 - Abb63 + Abb64
      AbbSumP469 = Abb64 + 1/3.D0*Abb18 - 1/3.D0*Abb63
      AbbSumP332 = Abb63 + AbbSumP85
      AbbSumP368 = AbbSumP66 - 1/2.D0*Abb24
      AbbSumP200 = Abb18 + Abb66 - 3/2.D0*Abb63
      AbbSumP286 = Abb18 - Abb63 + Abb66
      AbbSumP463 = Abb39 - Abb40 - 2*Abb69
      AbbSumP161 = AbbSumP276 - 1/2.D0*Abb69
      AbbSumP304 = Abb33 + AbbSumP221
      AbbSumP496 = 3*Abb39 - 3*Abb69 + Abb70
      AbbSumP400 = Abb70 + AbbSumP205
      AbbSumP471 = Abb70 + 1/3.D0*Abb39 - 1/3.D0*Abb69
      AbbSumP193 = Abb70 - AbbSumP205
      AbbSumP465 = Abb46 - Abb20 - 2*Abb62
      AbbSumP170 = AbbSumP280 - 1/2.D0*Abb62
      AbbSumP306 = Abb41 + AbbSumP224
      AbbSumP499 = 3*Abb46 - 3*Abb62 + Abb72
      AbbSumP407 = Abb72 + AbbSumP209
      AbbSumP473 = Abb72 + 1/3.D0*Abb46 - 1/3.D0*Abb62
      AbbSumP195 = Abb72 - AbbSumP209
      AbbSumP174 = -Abb82 + AbbSumP281
      AbbSumP307 = Abb45 + AbbSumP226
      AbbSumP408 = Abb32 - Abb82 + Abb74
      AbbSumP474 = Abb74 + 1/3.D0*Abb32 - 1/3.D0*Abb82
      AbbSumP338 = Abb82 + AbbSumP96
      AbbSumP371 = AbbSumP74 - 1/2.D0*Abb45
      AbbSumP210 = Abb32 + Abb76 - 3/2.D0*Abb82
      AbbSumP291 = Abb32 - Abb82 + Abb76
      AbbSumP4 = -Abb24 + AbbSumP349
      AbbSumP16 = -Abb33 + AbbSumP356
      AbbSumP27 = -Abb41 + AbbSumP361
      AbbSumP34 = -Abb45 + AbbSumP363
      AbbSumP452 = Abb64 + AbbSumP200
      AbbSumP59 = -AbbSumP288 + AbbSumP70
      AbbSumP62 = -AbbSumP290 + AbbSumP73
      AbbSumP456 = Abb74 + AbbSumP210
      AbbSumP236 = AbbSumP107 - 1/2.D0*Abb59 + 1/2.D0*Abb18
      AbbSumP319 = Abb80 - AbbSumP286
      AbbSumP237 = AbbSumP110 - 1/2.D0*Abb53 + 1/2.D0*Abb32
      AbbSumP320 = Abb73 - AbbSumP291
      AbbSumP437 = Abb29 + Abb65 - Abb66 - 1/2.D0*Abb24 + 3/2.D0*Abb63
      AbbSumP56 = Abb65 - Abb66 + AbbSumP332
      AbbSumP422 = 2*Abb39 - Abb40 - 3*Abb69 + AbbSumP221
      AbbSumP312 = AbbSumP70 - 1/2.D0*Abb50 - 1/2.D0*Abb39 + 1/2.D0*Abb69
      AbbSumP429 = 2*Abb46 - Abb20 - 3*Abb62 + AbbSumP224
      AbbSumP313 = AbbSumP73 - 1/2.D0*Abb55 - 1/2.D0*Abb46 + 1/2.D0*Abb62
      AbbSumP444 = Abb56 + Abb75 - Abb76 - 1/2.D0*Abb45 + 3/2.D0*Abb82
      AbbSumP63 = Abb75 - Abb76 + AbbSumP338
      AbbSumP385 = Abb79 - AbbSumP368 - 1/2.D0*Abb59
      AbbSumP487 = -Abb64 - Abb66 + Abb80 - 3/2.D0*Abb18 + 3/2.D0*Abb63
      AbbSumP386 = Abb81 - AbbSumP371 - 1/2.D0*Abb53
      AbbSumP488 = -Abb74 - Abb76 + Abb73 - 3/2.D0*Abb32 + 3/2.D0*Abb82
      AbbSumP191 = Abb65 - Abb66 + AbbSumP85 - 1/2.D0*Abb24 + 3/2.D0*Abb63
      AbbSumP447 = -Abb65 + AbbSumP452 + 1/2.D0*Abb24
      AbbSumP196 = Abb75 - Abb76 + AbbSumP96 - 1/2.D0*Abb45 + 3/2.D0*Abb82
      AbbSumP450 = -Abb75 + AbbSumP456 + 1/2.D0*Abb45
      AbbSumP311 = AbbSumP107 - AbbSumP56
      AbbSumP314 = AbbSumP110 - AbbSumP63
      AbbSumP414 = -2*Abb18 + Abb29 + 4*Abb63 + 2*Abb65 - 
     -     2*Abb66 - AbbSumP302
      AbbSumP432 = -2*Abb32 + Abb56 + 4*Abb82 + 2*Abb75 - 
     -     2*Abb76 - AbbSumP307
      AbbSumP505 = -Abb65 + Abb66 + AbbSumP107 - 1/2.D0*Abb59 + 
     -     3/2.D0*Abb18 + 1/2.D0*Abb24 - 3/2.D0*Abb63
      AbbSumP510 = -Abb75 + Abb76 + AbbSumP110 - 1/2.D0*Abb53 + 
     -     3/2.D0*Abb32 + 1/2.D0*Abb45 - 3/2.D0*Abb82
      AbbSumP345 = Abb66 + AbbSumP107 - AbbSumP368 - 1/2.D0*Abb59 + 
     -     3/2.D0*Abb18 - 3/2.D0*Abb63
      AbbSumP491 = 2*Abb18 - 2*Abb63 - Abb65 + Abb66 + 
     -     AbbSumP107 - 1/2.D0*Abb59 + 1/2.D0*Abb24 + 1/2.D0*Abb64
      AbbSumP476 = Abb64 + AbbSumP505
      AbbSumP459 = Abb24 - Abb29 - 3*Abb63 - 2*Abb65 + 
     -     2*Abb66 + AbbSumP107 - 1/2.D0*Abb59 + 3/2.D0*Abb18
      AbbSumP346 = Abb76 + AbbSumP110 - AbbSumP371 - 1/2.D0*Abb53 + 
     -     3/2.D0*Abb32 - 3/2.D0*Abb82
      AbbSumP501 = 2*Abb32 - 2*Abb82 - Abb75 + Abb76 + 
     -     AbbSumP110 - 1/2.D0*Abb53 + 1/2.D0*Abb45 + 1/2.D0*Abb74
      AbbSumP477 = Abb74 + AbbSumP510
      AbbSumP466 = Abb45 - Abb56 - 3*Abb82 - 2*Abb75 + 
     -     2*Abb76 + AbbSumP110 - 1/2.D0*Abb53 + 3/2.D0*Abb32

** Abbreviation sums for boxes attached to the lower quark line
      AbbSumL30 = Abb55 - Abb46
      AbbSumL23 = Abb59 - Abb18
      AbbSumL70 = Abb1 - 4*Abb2
      AbbSumL75 = Abb1 - 1/2.D0*Abb2
      AbbSumL95 = Abb1 + 2*Abb2
      AbbSumL100 = Abb1 + 3*Abb2
      AbbSumL33 = Abb53 - Abb32
      AbbSumL27 = Abb50 - Abb39
      AbbSumL1 = Abb1 + Abb3
      AbbSumL46 = Abb41 - Abb62
      AbbSumL43 = Abb24 - Abb63
      AbbSumL47 = Abb45 - Abb82
      AbbSumL44 = Abb33 - Abb69
      AbbSumL2 = Abb2 + Abb4
      AbbSumL76 = -2*Abb3 + Abb4
      AbbSumL72 = Abb4 - 1/4.D0*Abb3
      AbbSumL101 = Abb4 + 1/3.D0*Abb3
      AbbSumL97 = Abb4 + 1/2.D0*Abb3
      AbbSumL4 = AbbSumL23 + AbbSumL30
      AbbSumL7 = AbbSumL27 + AbbSumL33
      AbbSumL55 = AbbSumL30 + AbbSumL46
      AbbSumL53 = AbbSumL23 + AbbSumL43
      AbbSumL12 = AbbSumL43 + AbbSumL46
      AbbSumL56 = AbbSumL33 + AbbSumL47
      AbbSumL87 = AbbSumL33 - 1/3.D0*Abb45 + 1/3.D0*Abb82
      AbbSumL116 = AbbSumL33 - 2/3.D0*Abb45 + 2/3.D0*Abb82
      AbbSumL74 = AbbSumL33 - AbbSumL47
      AbbSumL54 = AbbSumL27 + AbbSumL44
      AbbSumL13 = AbbSumL44 + AbbSumL47
      AbbSumL86 = AbbSumL27 - 1/3.D0*Abb33 + 1/3.D0*Abb69
      AbbSumL114 = AbbSumL27 - 2/3.D0*Abb33 + 2/3.D0*Abb69
      AbbSumL71 = AbbSumL27 - AbbSumL44
      AbbSumL66 = -2*Abb1 - 2*Abb3 + AbbSumL2
      AbbSumL63 = AbbSumL2 - 1/4.D0*Abb1 - 1/4.D0*Abb3
      AbbSumL92 = AbbSumL2 + 1/3.D0*Abb1 + 1/3.D0*Abb3
      AbbSumL90 = AbbSumL2 + 1/2.D0*Abb1 + 1/2.D0*Abb3
      AbbSumL16 = AbbSumL12 + AbbSumL4
      AbbSumL17 = AbbSumL13 + AbbSumL7
      AbbSumL79 = AbbSumL7 - 1/3.D0*Abb45 + 1/3.D0*Abb82 - 
     -     1/3.D0*Abb33 + 1/3.D0*Abb69
      AbbSumL108 = AbbSumL7 - 2/3.D0*Abb45 + 2/3.D0*Abb82 - 
     -     2/3.D0*Abb33 + 2/3.D0*Abb69
      AbbSumL65 = -AbbSumL13 + AbbSumL7

** Abbreviation sums for boxes attached to the upper quark line
      AbbSumU3 = Abb1 + Abb2
      AbbSumU37 = -Abb62 + Abb41
      AbbSumU123 = Abb41 - 1/2.D0*Abb62
      AbbSumU41 = Abb41 + 1/2.D0*Abb62
      AbbSumU72 = 2*Abb62 + Abb41
      AbbSumU50 = Abb45 + Abb41
      AbbSumU79 = Abb33 + Abb24
      AbbSumU118 = Abb33 - Abb69
      AbbSumU122 = Abb33 - 1/2.D0*Abb69
      AbbSumU11 = Abb1 - Abb3
      AbbSumU130 = Abb1 + 1/4.D0*Abb3
      AbbSumU99 = Abb1 + 1/3.D0*Abb3
      AbbSumU63 = Abb1 + 2*Abb3
      AbbSumU27 = Abb1 + 7/3.D0*Abb3
      AbbSumU105 = Abb1 + 3*Abb3
      AbbSumU133 = Abb1 + 4*Abb3
      AbbSumU117 = Abb24 - Abb63
      AbbSumU120 = Abb24 - 1/2.D0*Abb63
      AbbSumU12 = -Abb2 + Abb4
      AbbSumU135 = Abb4 + 1/4.D0*Abb2
      AbbSumU107 = Abb4 + 1/3.D0*Abb2
      AbbSumU29 = Abb4 + 3/7.D0*Abb2
      AbbSumU64 = Abb4 + 1/2.D0*Abb2
      AbbSumU100 = 3*Abb2 + Abb4
      AbbSumU131 = 4*Abb2 + Abb4
      AbbSumU1 = Abb3 + Abb4
      AbbSumU39 = Abb45 - Abb82
      AbbSumU124 = Abb45 - 1/2.D0*Abb82
      AbbSumU43 = Abb45 + 1/2.D0*Abb82
      AbbSumU4 = Abb62 + Abb82
      AbbSumU74 = Abb45 + 2*Abb82
      AbbSumU86 = -Abb69 - Abb63 + AbbSumU79
      AbbSumU88 = AbbSumU79 - 1/2.D0*Abb69 - 1/2.D0*Abb63
      AbbSumU5 = AbbSumU1 - AbbSumU3
      AbbSumU96 = Abb3 + AbbSumU135 + 1/4.D0*Abb1
      AbbSumU81 = Abb3 + AbbSumU107 + 1/3.D0*Abb1
      AbbSumU16 = Abb3 + AbbSumU29 + 3/7.D0*Abb1
      AbbSumU55 = Abb3 + AbbSumU64 + 1/2.D0*Abb1
      AbbSumU77 = 3*Abb1 + Abb3 + AbbSumU100
      AbbSumU94 = 4*Abb1 + Abb3 + AbbSumU131
      AbbSumU20 = AbbSumU37 + AbbSumU39
      AbbSumU89 = AbbSumU123 + AbbSumU124
      AbbSumU22 = AbbSumU41 + AbbSumU43
      AbbSumU60 = AbbSumU72 + AbbSumU74

** Abbreviation sums for the corrections to the qqV vertex
      AbbSumQ35 = Abb1 - Abb2
      AbbSumQ60 = Abb1 - 1/3.D0*Abb2
      AbbSumQ9 = Abb1 + 1/3.D0*Abb3
      AbbSumQ1 = Abb1 + 1/2.D0*Abb3
      AbbSumQ25 = Abb1 + 5/3.D0*Abb3
      AbbSumQ5 = Abb1 + 2*Abb3
      AbbSumQ15 = Abb1 + 3*Abb3
      AbbSumQ17 = Abb4 + 1/3.D0*Abb2
      AbbSumQ7 = Abb4 + 1/2.D0*Abb2
      AbbSumQ26 = Abb4 + 3/5.D0*Abb2
      AbbSumQ3 = 2*Abb2 + Abb4
      AbbSumQ12 = 3*Abb2 + Abb4
      AbbSumQ61 = -3*Abb3 + Abb4
      AbbSumQ40 = -Abb3 + Abb4


      end


********************************************************************************
