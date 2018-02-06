********************************************************************************
********************************************************************************
***                                                                          ***
*** mat_elNC+.F                                                              ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This subroutine (and those accompanying it) calculate the squared matrix ***
*** elements for the boxes and pentagons for the charged currents. Note that ***
*** we're using FormCalc notation.  Much of this is purloined from           ***
*** FormCalc's SquaredME.F.                                                  ***
***                                                                          ***
********************************************************************************
********************************************************************************

      subroutine squaringCC(p,mass,qk,CC,SInv,TInv,SqMatEl,matFAC,debug)

      implicit none


#include "global.inc"

** Particle label and do-loop variable
**        i  =  1         upper incoming
**           =  2         upper outgoing
**           =  3         lower incoming
**           =  4         lower outgoing
**           =  5         higgs outgoing
      integer i, Gen6


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
** NOTE that, by default, the qqV loops are not calculated - instead we use an
** effective qqV coupling as it's faster
      double complex CloopNlp(1), CloopNup(1), CloopPPp(1)
      double complex CloopPBl(1), CloopPBu(1)
      double complex CloopNlnop(1), CloopNunop(1), CloopPPnop(1) 
      double complex CloopCT(1), CloopUPp(1), CloopUPn(1)
      double complex CloopDOWNn(1), CloopDOWNp(1)

                                                    
** The squared matrix element. The first entry is current, the second is tree 
** (1) or loop (2) level
      double precision SqMatEl(6,2), SubMatEl


** 4-momenta p(mu,i) of the external quarks.  
      double precision p(0:3,5)


** Subtraction term from real corrections:  see notebook 7, page 14
      double precision Dsub


** Charged current parameter: see notebook 7, page 3.  Cur = 4 + CC
**      CC = 1     u + s -> d + H + c
**         = 2     d + c -> u + H + c       
      integer CC, cur

** qk is the quark/antiquark pattern label. It is dependent on fsign: 
**               qk = 1          quark-quark
**                  = 2          antiquark-quark
**                  = 3          quark-antiquark
**                  = 4          antiquark-antiquark
      integer qk


** Dummy parameter
      integer dummy

** switch deciding whether to run qqV matrix elements: false by default
** NOTE that, if you switch this to .true. you should also set qqVmat to true
** in mat_elNC+.F, and comment out the call to qqV_eff in 
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


** These declare the abbreviations, the integral id's and the parameters
#include "abbr_CC.h"
#include "integrals_CC.h"
#include "mssm.inc"

** Needed for looptools
#include "looptools.h"


      if (debug .gt. 0) then
         write(*,*)'CC =', CC
      end if


      cur = CC + 4


** Sorting out the SUSY factors
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

      if (debug .ne. 0 .and. CC .eq. 1 .and. qk .eq. 1) then
         write(*,*)'S, T, T14, U, T24, S34 =', S, T, T14, U, T24, S34
*         stop
      end if


      MatSUN = 9D0


** Setting the current-specific parameters
      do i = 1, 5
         mass2(i) = mass(i)**2
         q(i) = chargeCur(cur,i)
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


** Debugging output
      if (debug .ne. 0) then
         do i = 1, 4
            write(*,*)'i =', i
            write(*,*)'mass(i) =', mass(i)
            write(*,*)'q(i) =', q(i)
            write(*,*)'gL(i), gR(i) =', gL(i), gR(i)
            write(*,*)'p(0,i), p(1,i) =', p(0,i), p(1,i)
            write(*,*)'p(2,i), p(3,i) =', p(2,i), p(3,i)
            write(*,*)'   '            
         end do
      end if 


** The real corrections
      Dsub = (0D0,0D0) 
      call realCor(mass,qk,cur,p,q,Dsub,debug)


** The integral values are calculated: suffix L = lower box, U = upper box,
** P = pentagon
** C integrals ...
      iintL1 = Cget(T24,T,mass2(4),MZ2,mass2(3),mass2(4))
      iintL3 = Cget(mass2(4),T24,T,0.D0,mass2(4),MW2)
      iintL4 = Cget(T,T24,mass2(4),0.D0,MW2,mass2(4))
      iintL5 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,T,0.D0,MW2,MW2)
      iintL6 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,T,0.D0,MW2,MW2)
      iintL11 = Cget(mass2(4),T24,T,MZ2,mass2(4),MW2)
      iintL12 = Cget(MH2,T,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,MW2,MW2,MZ2)
      iintL13 = Cget(T,mass2(4),T24,MW2,MZ2,mass2(4))
      iintL14 = Cget(T,T24,mass2(4),MW2,MZ2,mass2(3))
      iintL15 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,T,MW2,MZ2,MZ2)
      iintL16 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,T,MH2,MW2,MZ2,MW2)


      iintU1 = Cget(T,MH2,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,0.D0,MW2,MW2)
      iintU2 = Cget(T,MH2,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,0.D0,MW2,MW2)
      iintU3 = iintL6
      iintU8 = iintL12
      iintU9 = Cget(MH2,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,T,MW2,MW2,MZ2)
      iintU10 = Cget(T,MH2,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,MW2,MZ2,MZ2)
      iintU11 = Cget(T,mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,MW2,MZ2,MW2)
      iintU12 = iintL15
      iintU19 = Cget(T14,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(2),MZ2,mass2(1),mass2(2))

** D integrals ...
      iintL2 = Dget(MH2,T24,T,
     -     mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -     T14 - T24,mass2(3),mass2(4),MZ2,MZ2,mass2(3),mass2(4))
      iintL7 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,mass2(3),T24,T,mass2(4),MH2,
     -     0.D0,MW2,mass2(4),MW2)
      iintL8 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,mass2(4),
     -     mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -     T14 - T24,T,mass2(3),MH2,0.D0,MW2,mass2(3),MW2)
      iintL9 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,
     -     mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -     T14 - T24,mass2(3),T,mass2(4),0.D0,MW2,MW2,mass2(3))
      iintL10 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,T24,mass2(4),T,mass2(3),
     -     0.D0,MW2,MW2,mass2(4))
      iintL17 = Dget(MH2,T,mass2(3),mass2(4),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(1) + mass2(2) + mass2(3) + mass2(4) + 
     -     2*MH2 - S34 - T14 - T24,MW2,MW2,MZ2,mass2(3))
      iintL18 = Dget(MH2,T,mass2(4),mass2(3),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,T24,MW2,MW2,MZ2,mass2(4))
      iintL19 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,mass2(3),
     -     mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -     T14 - T24,MH2,mass2(4),T,MW2,MZ2,mass2(3),MW2)
      iintL20 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,mass2(4),T24,MH2,mass2(3),T,
     -     MW2,MZ2,mass2(4),MW2)
      iintL21 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,
     -     mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -     T14 - T24,mass2(3),T,mass2(4),MW2,MZ2,MZ2,mass2(4))
      iintL22 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,T24,mass2(4),T,mass2(3),
     -     MW2,MZ2,MZ2,mass2(3))

      iintU4 = Dget(T,mass2(1),T14,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(2),MH2,0.D0,MW2,mass2(2),MW2)
      iintU5 = Dget(T,mass2(2),S34,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(1),MH2,0.D0,MW2,mass2(1),MW2)
      iintU6 = Dget(T,MH2,S34,mass2(1),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(2),0.D0,MW2,MW2,mass2(1))
      iintU7 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(1),mass2(2),T,T14,
     -     0.D0,MW2,MW2,mass2(2))
      iintU13 = Dget(MH2,T,mass2(2),T14,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(1),MW2,MW2,MZ2,mass2(2))
      iintU14 = Dget(MH2,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(1),mass2(2),T,S34,MW2,MW2,MZ2,mass2(1))
      iintU15 = Dget(T,mass2(1),S34,MH2,mass2(2),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,MW2,MZ2,mass2(1),MW2)
      iintU16 = Dget(T,mass2(2),T14,MH2,mass2(1),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,MW2,MZ2,mass2(2),MW2)
      iintU17 = Dget(T,MH2,S34,mass2(1),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(2),MW2,MZ2,MZ2,mass2(2))
      iintU18 = Dget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(1),mass2(2),T,T14,
     -     MW2,MZ2,MZ2,mass2(1))
      iintU20 = Dget(MH2,T14,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,S34,mass2(1),mass2(2),MZ2,MZ2,mass2(1),mass2(2))


      if (qk .eq. 1 .or. qk .eq. 4) then
         iintP1 = Dget(T,MH2,mass2(3),mass2(4),
     -        mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -        T24 - U,T24,0.D0,MW2,MW2,mass2(4))
         iintP2 = iintL9
         iintP7 = iintL17
         iintP8 = Dget(MH2,
     -        mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -        T24 - U,mass2(4),T24,T,mass2(3),MW2,MW2,MZ2,mass2(4))
         iintP9 = Dget(T,MH2,mass2(3),mass2(4),
     -        mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -        T24 - U,T24,MW2,MZ2,MZ2,mass2(3))
         iintP10 = iintL21
      else if (qk .eq. 2 .or. qk .eq. 3) then
         iintP1 = Dget(MH2,mass2(2),S,mass2(4),S34,
     -        mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -        T14 - T24,MW2,MW2,mass2(1),mass2(3))
         iintP2 = iintL10
         iintP7 = Dget(MH2,mass2(2),S,mass2(4),S34,
     -        mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -        T14 - T24,MZ2,MZ2,mass2(2),mass2(4))
         iintP8 = iintL18
         iintP9 = iintL22
      end if
      
** E integrals ...
      iintP3 = Eget(T,MH2,S34,
     -     2*mass2(1) + mass2(2) + mass2(3) + mass2(4) + MH2 - S - 
     -     T - T14,mass2(4),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(2),mass2(3),mass2(1),T24,0.D0,MW2,MW2,
     -     mass2(1),mass2(4))
      iintP4 = Eget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(1),
     -     mass2(1) + 2*mass2(2) + mass2(3) + mass2(4) + MH2 - S34 - 
     -     T - U,mass2(4),T,T14,T24,mass2(2),mass2(3),0.D0,MW2,MW2,
     -     mass2(2),mass2(4))
      iintP5 = Eget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(1),U,mass2(3),T,
     -     T14,mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - 
     -     S34 - T14 - T24,mass2(2),mass2(4),0.D0,MW2,MW2,mass2(2),
     -     mass2(3))
      iintP6 = Eget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(2),S,mass2(3),T,
     -     S34,mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - 
     -     S34 - T14 - T24,mass2(1),mass2(4),0.D0,MW2,MW2,mass2(1),
     -     mass2(3))
      iintP11 = Eget(MH2,T,mass2(1),S,mass2(4),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(2),mass2(3),S34,
     -     mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -     T14 - T24,MW2,MW2,MZ2,mass2(1),mass2(3))
      iintP12 = Eget(MH2,T,mass2(2),
     -     mass2(1) + 2*mass2(2) + mass2(3) + mass2(4) + MH2 - S34 - 
     -     T - U,mass2(3),mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,mass2(1),mass2(4),T14,T24,
     -     MW2,MW2,MZ2,mass2(2),mass2(4))
      iintP13 = Eget(MH2,T,mass2(2),U,mass2(4),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(1),mass2(3),T14,
     -     mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - S34 - 
     -     T14 - T24,MW2,MW2,MZ2,mass2(2),mass2(3))
      iintP14 = Eget(MH2,
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(1),
     -     2*mass2(1) + mass2(2) + mass2(3) + mass2(4) + MH2 - S - 
     -     T - T14,T24,T,S34,mass2(4),mass2(2),mass2(3),MW2,MW2,MZ2,
     -     mass2(1),mass2(4))
      iintP15 = Eget(T,MH2,S34,
     -     2*mass2(1) + mass2(2) + mass2(3) + mass2(4) + MH2 - S - 
     -     T - T14,mass2(4),
     -     mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -     T24 - U,mass2(2),mass2(3),mass2(1),T24,MW2,MZ2,MZ2,
     -     mass2(2),mass2(3))
      iintP16 = Eget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(1),
     -     mass2(1) + 2*mass2(2) + mass2(3) + mass2(4) + MH2 - S34 - 
     -     T - U,mass2(4),T,T14,T24,mass2(2),mass2(3),MW2,MZ2,MZ2,
     -     mass2(1),mass2(3))
      iintP17 = Eget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(1),U,mass2(3),T,
     -     T14,mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - 
     -     S34 - T14 - T24,mass2(2),mass2(4),MW2,MZ2,MZ2,mass2(1),
     -     mass2(4))
      iintP18 = Eget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - S - T24 - U,MH2,mass2(2),S,mass2(3),T,
     -     S34,mass2(1) + mass2(2) + mass2(3) + mass2(4) + 2*MH2 - 
     -     S34 - T14 - T24,mass2(1),mass2(4),MW2,MZ2,MZ2,mass2(2),
     -     mass2(4))

      if (qk .eq. 2 .or. qk .eq. 3) then
         iintP10 = iintP11
         iintP11 = iintP12
         iintP12 = iintP13
         iintP13 = iintP14
         iintP14 = iintP15
         iintP15 = iintP16
         iintP16 = iintP17
         iintP17 = iintP18
      end if


** qqV integrals: (these are called iff we're not using an effective coupling to
** take account of the corrections to the qqV vertex)
      if (qqVmat) then     
         iintQl1 = Bget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -        mass2(4) + MH2 - S - T24 - U,mass2(3),mass2(4))
         iintQl2 = Cget(mass2(3),
     -        mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -        T24 - U,mass2(4),0.D0,mass2(3),mass2(4))
         iintQl3 = Cget(mass2(3),
     -        mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - S - 
     -        T24 - U,mass2(4),MZ2,mass2(3),mass2(4))
         iintQl4 = Bget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -        mass2(4) + MH2 - S - T24 - U,0.D0,MW2)
         iintQl5 = Bget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -        mass2(4) + MH2 - S - T24 - U,MW2,MZ2)
         iintQl6 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -        mass2(4) + MH2 - S - T24 - U,mass2(3),mass2(4),0.D0,MW2,
     -        mass2(4))
         iintQl7 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -        mass2(4) + MH2 - S - T24 - U,mass2(4),mass2(3),MW2,MZ2,
     -        mass2(4))
         iintQl8 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -        mass2(4) + MH2 - S - T24 - U,mass2(3),mass2(4),MW2,MZ2,
     -        mass2(3))
         iintQl9 = Cget(mass2(1) + mass2(2) + 2*mass2(3) + 
     -        mass2(4) + MH2 - S - T24 - U,mass2(4),mass2(3),0.D0,MW2,
     -        mass2(3))
         
         iintQu1 = Bget(T,mass2(1),mass2(2))
         iintQu2 = Cget(mass2(1),T,mass2(2),0.D0,mass2(1),mass2(2))
         iintQu3 = Cget(mass2(1),T,mass2(2),MZ2,mass2(1),mass2(2))
         iintQu4 = Bget(T,0.D0,MW2)
         iintQu5 = Bget(T,MW2,MZ2)
         iintQu6 = Cget(T,mass2(1),mass2(2),MW2,MZ2,mass2(1))
         iintQu7 = Cget(T,mass2(2),mass2(1),0.D0,MW2,mass2(1))
         iintQu8 = Cget(T,mass2(1),mass2(2),0.D0,MW2,mass2(2))
         iintQu9 = Cget(T,mass2(2),mass2(1),MW2,MZ2,mass2(2))
         

         if (qk .eq. 2 .or. qk .eq. 4) then
            dummy = iintQu6
            iintQu6 = iintQu8
            iintQu8 = dummy
            
            dummy = iintQu7
            iintQu7 = iintQu9
            iintQu9 = dummy
         end if
         if (qk .ge. 3) then
            dummy = iintQl6
            iintQl6 = iintQl8
            iintQl8 = dummy
            
            dummy = iintQl7
            iintQl7 = iintQl9
            iintQl9 = dummy
         end if
      end if


** Clearing the values of the squared matrix element
      SubMatEl = 0D0
      SqMatEl(CC+4,1) = 0D0
      SqMatEl(CC+4,2) = 0D0


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
            call hel_qqCC 
         else if (qk .eq. 2) then
            call hel_qbarqCC
         else if (qk .eq. 3) then
            call hel_qqbarCC
         else if (qk .eq. 4) then
            call hel_qbarqbarCC
         end if


** Clearing the values of the tree and loop components
         Ctree(1) = (0D0,0D0)
         Cloop(1) = (0D0,0D0)
         
         CloopCT(1) = (0D0,0D0)
         CloopUPp(1) = (0D0,0D0)
         CloopUPn(1) = (0D0,0D0)
         CloopDOWNp(1) = (0D0,0D0)
         CloopDOWNn(1) = (0D0,0D0)

         CloopNlp(1) = (0D0,0D0)
         CloopNup(1) = (0D0,0D0)
         CloopPPp(1) = (0D0,0D0)
         CloopPBl(1) = (0D0,0D0)

         CloopNlnop(1) = (0D0,0D0)
         CloopNunop(1) = (0D0,0D0)
         CloopPPnop(1) = (0D0,0D0)
         CloopPBu(1) = (0D0,0D0)
                          
                          
** Tree level
         call treeCC(mass2,S,T,U,T24,Ctree)


** The pentagons
         if (qk .eq. 1) then
            call pentP_p_qq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPPp)
            call pentP_nop_qq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPPnop)
         else if (qk .eq. 2) then
            call pentP_p_qbarq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPPp)
            call pentP_nop_qbarq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -           gR,CloopPPnop)
         else if (qk .eq. 3) then
            call pentP_p_qqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPPp)
            call pentP_nop_qqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -           gR,CloopPPnop)
         else if (qk .eq. 4) then
            call pentP_p_qbarqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -           gR,CloopPPp)
            call pentP_nop_qbarqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,
     -           gL,gR,CloopPPnop)
         end if
         
** The boxes:
         if (qk .eq. 1 .or. qk .eq. 3) then
            call boxesNu_p13_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopNup)
            call boxesNu_nop13_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -           gR,CloopNunop)
            call boxesPu_13_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPBu)
         else if (qk .eq. 2 .or. qk .eq. 4) then
            call boxesNu_p24_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopNup)
            call boxesNu_nop24_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -           gR,CloopNunop)
            call boxesPu_24_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPBu)
         end if
         
         if (qk .eq. 1 .or. qk .eq. 2) then
            call boxesNl_p12_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopNlp)
            call boxesNl_nop12_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -           gR,CloopNlnop)
            call boxesPl_12_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPBl)
         else if (qk .eq. 3 .or. qk .eq. 4) then
            call boxesNl_p34_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopNlp)
            call boxesNl_nop34_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -           gR,CloopNlnop)
            call boxesPl_34_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -           CloopPBl)
         end if
                          
** The qqV vertex: (this is called iff we're not using an effective qqV coupling
         if (qqVmat) then
            call qqVCT(qk,CC,mass2,S,T,T24,T14,S34,U,CloopCT)
            
            if (qk .eq. 1 .or. qk .eq. 3) then
               call qqVUP_n13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopUPn)
               call qqVUP_p13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopUPp)
            else if (qk .eq. 2 .or. qk .eq. 4) then
               call qqVUP_n24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopUPn)
               call qqVUP_p24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopUPp)
            end if
            
            if (qk .le. 2) then
               call qqVDOWN_n12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopDOWNn)
               call qqVDOWN_p12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopDOWNp)
            else if (qk .ge. 3) then
               call qqVDOWN_n34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopDOWNn)
               call qqVDOWN_p34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -              CloopDOWNp)
            end if
         end if


**Set individual components to zero, for debugging purposes
*         CloopCT(1) = (0D0,0D0)
*         CloopUPn(1) = (0D0,0D0)
*         CloopUPp(1) = (0D0,0D0)
*         CloopDOWNn(1) = (0D0,0D0)
*         CloopDOWNp(1) = (0D0,0D0)
         
*         CloopNlp(1) = (0D0,0D0)
*         CloopNup(1) = (0D0,0D0)
*         CloopPPp(1) = (0D0,0D0)
*         CloopPBl(1) = (0D0,0D0)
         
*         CloopNlnop(1) = (0D0,0D0)
*         CloopNunop(1) = (0D0,0D0)
*         CloopPPnop(1) = (0D0,0D0)
*         CloopPBu(1) = (0D0,0D0)


** Correcting the signs of each diagram set: for the qqV vertex, then boxes
         if (CC .eq. 2) then
            CloopUPn(1) = -CloopUPn(1)
            CloopDOWNn(1) = -CloopDOWNn(1)
            
            CloopNlp(1) = -CloopNlp(1)
            CloopNup(1) = -CloopNup(1)
            CloopNlnop(1) = -CloopNlnop(1)
            CloopNunop(1) = -CloopNunop(1)
         end if
                              
         if (qk .eq. 2) then
            Ctree(1) = -Ctree(1)
            
            CloopCT(1) = -CloopCT(1)
            CloopDOWNn(1) = -CloopDOWNn(1)
            CloopDOWNp(1) = -CloopDOWNp(1)
            CloopUPp(1) = -CloopUPp(1)
            CloopUPn(1) = -CloopUPn(1)
            
            CloopNlp(1) = -CloopNlp(1)
            CloopNlnop(1) = -CloopNlnop(1)
            CloopPBl(1) = -CloopPBl(1)
** Matrix element has been pre-multiplied by -1
*            CloopPBu(1) = -CloopPBu(1)
         end if
                              
         if (qk .eq. 3) then
            Ctree(1) = -Ctree(1)
            
            CloopCT(1) = -CloopCT(1)
            CloopUPp(1) = -CloopUPp(1)
            CloopUPn(1) = -CloopUPn(1)
            CloopDOWNn(1) = -CloopDOWNn(1)
            CloopDOWNp(1) = -CloopDOWNp(1)
            
            CloopNup(1) = -CloopNup(1)
            CloopNunop(1) = -CloopNunop(1)
** Matrix element has been pre-multiplied by -1
*            CloopPBl(1) = -CloopPBl(1)
            CloopPBu(1) = -CloopPBu(1)
         end if

         if (qk .eq. 4) then
            CloopNlp(1) = -CloopNlp(1)
            CloopNlnop(1) = -CloopNlnop(1)
            CloopNup(1) = -CloopNup(1)
            CloopNunop(1) = -CloopNunop(1)
** Added because of pre-multiplication of matrix element by -1
            CloopPBl(1) = -CloopPBl(1)
            CloopPBu(1) = -CloopPBu(1)
         end if
                          

** The total loop component 
         Cloop(1) = CloopNlp(1) + CloopNup(1) + CloopNlnop(1) + 
     -        CloopNunop(1) + CloopPPp(1) + CloopPPnop(1) + 
     -        CloopPBl(1) + CloopPBu(1) + CloopCT(1) + CloopUPp(1) + 
     -        CloopUPn(1) + CloopDOWNp(1) + CloopDOWNn(1)


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


** The squared matrix elements
         SqMatEl(CC+4,1) = SqMatEl(CC+4,1) + sumup(MatSUN,Ctree,Ctree)  
                              
         SqMatEl(CC+4,2) = SqMatEl(CC+4,2) + 2*sumup(MatSUN,Cloop,Ctree)  


** Debugging messages 
         if (debug .ne. 0) then
            write(*,*)'Hel(i) =', Hel
            
            if (debug .gt. 1) then
               write(*,*)'Ctree(1) =', Ctree(1)
               
               
               write(*,*)'CloopCT =', CloopCT(1)  
               write(*,*)'CloopUPp =', CloopCT(1) + CloopUPp(1) 
               write(*,*)'CloopUPn =', CloopCT(1) + CloopUPp(1) + 
     -              CloopUPn(1)
               write(*,*)'CloopDOWNp =', CloopCT(1) + CloopUPp(1) +  
     -              CloopUPn(1) + CloopDOWNp(1)
               write(*,*)'CloopDOWNn =', CloopCT(1) + CloopUPp(1) + 
     -              CloopUPn(1) + CloopDOWNp(1) + CloopDOWNn(1)
               
               write(*,*)'CloopNlp =', CloopNlp(1)
               write(*,*)'CloopNlnop =', CloopNlnop(1) + CloopNlp(1)
               write(*,*)'CloopNup =', CloopNlnop(1) + CloopNlp(1) + 
     -              CloopNup(1)
               write(*,*)'CloopNunop =', CloopNlnop(1) + CloopNlp(1) + 
     -              CloopNup(1) + CloopNunop(1)
               write(*,*)'CloopPBu =', CloopNlnop(1) + CloopNlp(1) + 
     -              CloopNup(1) + CloopNunop(1) + CloopPBu(1)
               write(*,*)'CloopPBl =', CloopNlnop(1) + CloopNlp(1) + 
     -              CloopNup(1) + CloopNunop(1) + CloopPBu(1) +
     -              CloopPBl(1)
               write(*,*)'CloopPPp =',  CloopNlnop(1) + CloopNlp(1) + 
     -              CloopNup(1) + CloopNunop(1) + CloopPBu(1) +
     -              CloopPBl(1) + CloopPPp(1)
               write(*,*)'CloopPPnop =',  CloopNlnop(1) + CloopNlp(1) + 
     -              CloopNup(1) + CloopNunop(1) + CloopPBu(1) + 
     -              CloopPBl(1) + CloopPPp(1) + CloopPPnop(1)
                                    
            end if
                           
            write(*,*)'SqMatEl(1) =', SqMatEl(CC+4,1)
            write(*,*)'SqMatEl(2) =', SqMatEl(CC+4,2)
            write(*,*)'Dsub =', Dsub
            
         end if 
                              

** This will stop the code and output debugging messages if we get silly 
** answers for the loop contribution
         if (SqMatEl(CC+4,2) .eq. 1000D0) then
            write(*,*)'oops CC, qk =', CC, qk
            write(*,*)'Ctree(1) =', Ctree(1)
            
            write(*,*)'CloopCT =', CloopCT(1)  
            write(*,*)'CloopUPp =', CloopCT(1) + CloopUPp(1) 
            write(*,*)'CloopUPn =', CloopCT(1) + CloopUPp(1) + 
     -           CloopUPn(1)
            write(*,*)'CloopDOWNp =', CloopCT(1) + CloopUPp(1) + 
     -           CloopUPn(1) + CloopDOWNp(1)
            write(*,*)'CloopDOWNn =', CloopCT(1) +  CloopUPp(1) +  
     -           CloopUPn(1) + CloopDOWNp(1) + CloopDOWNn(1)
            
            write(*,*)'CloopNlp =', CloopNlp(1)
            write(*,*)'CloopNlnop =', CloopNlnop(1) + CloopNlp(1)
            write(*,*)'CloopNup =', CloopNlnop(1) + CloopNlp(1) + 
     -           CloopNup(1)
            write(*,*)'CloopNunop =', CloopNlnop(1) + CloopNlp(1) + 
     -           CloopNup(1) + CloopNunop(1)
            write(*,*)'CloopPBu =', CloopNlnop(1) + CloopNlp(1) + 
     -           CloopNup(1) + CloopNunop(1) + CloopPBu(1)
            write(*,*)'CloopPBl =', CloopNlnop(1) + CloopNlp(1) + 
     -           CloopNup(1) + CloopNunop(1) + CloopPBu(1) +
     -           CloopPBl(1)
            write(*,*)'CloopPPp =',  CloopNlnop(1) + CloopNlp(1) + 
     -           CloopNup(1) + CloopNunop(1) + CloopPBu(1) +
     -           CloopPBl(1) + CloopPPp(1)
            write(*,*)'CloopPPnop =',  CloopNlnop(1) + CloopNlp(1) + 
     -           CloopNup(1) + CloopNunop(1) + CloopPBu(1) +
     -           CloopPBl(1) + CloopPPp(1) + CloopPPnop(1)
            
            write(*,*)'Cloop =', Cloop(1)
            
            write(*,*)'SqMatEl(1) =', SqMatEl(CC+4,1)
            write(*,*)'SqMatEl(2) =', SqMatEl(CC+4,2)
            write(*,*)'Dsub =', Dsub
            
            write(*,*)'S, T, T14, U, T24, S34 =', S, T, T14, U, T24, S34
            write(*,*)'Abb1 =', Abb1
            stop
         end if
      end do


** Debugging messages 
      if (debug .ne. 0) then
         write(*,*)'   '
         write(*,*)'CC, qk =', CC, qk
         write(*,*)'Dsub =', Dsub
         write(*,*)'SqMatEl(CC+4,1) =', SqMatEl(CC+4,1)
         write(*,*)'SqMatEl(CC+4,2) =', SqMatEl(CC+4,2)         
      end if
 

** Adding on the subtraction term to deal with IR and collinear divergences
      SqMatEl(CC+4,2) = SqMatEl(CC+4,2) + (SubMatEl*Dsub)


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the tree part of the calculation

      subroutine treeCC(mass2,S,T,U,T24,Ctree)

      implicit none


** Momentum invariants
      double precision S, T, U, T24

** Tree level component
      double complex Ctree(1)

** Mass^2 of external particles
      double precision mass2(5)


** Parameter and abbreviation declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"


      Ctree(1) = Ctree(1) - (4*AlfaQED*MW*Pi)/SW2*(Abb1*EL)/
     -     (SW*(-MW2 + T)*(mass2(4) + mass2(2) + MH2 + 2*mass2(3) + 
     -     mass2(1) - MW2 - S - T24 - U))


      end


*******************************************************************************
*******************************************************************************

***************************************************************************
** This set of subroutines takes care of the calculation of the vector   **
** and spinor products.                                                  **
***************************************************************************

*******************************************************************************
*******************************************************************************

*** Subroutine that calculates the spinor products for q + q -> q + H + q
*** (qk = 1)

      subroutine hel_qqCC

      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel


** Vector products
      double complex F4, F6, F17, F19, F8, F5, F18, F11, F2, F12, F9
      double complex F15, F20, F22, F7, F14, F30, F29, F1, F3, F10
      double complex F28, F33, F21, F26, F16, F27, F34, F31, F32
      double complex F13, F23, F24, F25


** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS


** Contains the declarations of the abbreviations
#include "abbr_CC.h"


      F2 = vSeS(DottedSpinor(3,1,7),DottedSpinor(5,1,7))
      F1 = vSeS(Spinor(1,1,7),Spinor(2,1,7))

      F6 = vSxS(DottedSpinor(3,1,7),Spinor(2,1,6))
      F19 = vSxS(DottedSpinor(3,1,6),Spinor(2,1,7))
      F8 = vSxS(DottedSpinor(5,1,7),Spinor(1,1,6))
      F18 = vSxS(DottedSpinor(5,1,6),Spinor(1,1,7))
      F12 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     Spinor(1,1,7)))
      F15 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(1,1,7)))
      F20 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(2,1,7)))
      F22 = vSeS(DottedSpinor(3,1,6),DottedSpinor(5,1,6))
      F14 = vSxS(DottedSpinor(5,1,7),vVxS(k(1),
     -     Spinor(2,1,7)))
      F30 = vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(1,1,7)))
      F29 = vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(2,1,7)))
      F3 = vSeS(Spinor(1,1,6),Spinor(2,1,6))
      F10 = vSxS(DottedSpinor(3,1,7),vVeS(k(1),
     -     DottedSpinor(5,1,6)))
      F28 = vSxS(DottedSpinor(3,1,7),vVeS(k(2),
     -     DottedSpinor(5,1,6)))
      F33 = vSxS(DottedSpinor(3,1,7),vVeS(k(3),
     -     DottedSpinor(5,1,6)))
      F21 = vSxS(DottedSpinor(3,1,7),vVeS(k(4),
     -     DottedSpinor(5,1,6)))
      F26 = vSxS(DottedSpinor(3,1,6),vBeS(k(1),
     -     DottedSpinor(5,1,7)))
      F16 = vSxS(DottedSpinor(3,1,6),vBeS(k(2),
     -     DottedSpinor(5,1,7)))
      F27 = vSxS(DottedSpinor(3,1,6),vBeS(k(3),
     -     DottedSpinor(5,1,7)))
      F34 = vSxS(DottedSpinor(3,1,6),vBeS(k(4),
     -     DottedSpinor(5,1,7)))
      F31 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     vBeS(k(1),DottedSpinor(5,1,7))))
      F32 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     vBeS(k(3),DottedSpinor(5,1,7))))
      F13 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     vBeS(k(4),DottedSpinor(5,1,7))))
      F23 = vSxS(DottedSpinor(3,1,7),vVxS(k(3),
     -     vBeS(k(1),DottedSpinor(5,1,7))))
      F24 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     vBeS(k(1),DottedSpinor(5,1,7))))
      F25 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     vBeS(k(3),DottedSpinor(5,1,7))))


      Abb1 = F1*F2
      Abb8 = F1*F13
      Abb9 = F12*F14
      Abb10 = F14*F15
      Abb23 = F1*F23
      Abb24 = F1*F24
      Abb25 = F1*F25
      Abb29 = F12*F29
      Abb32 = F1*F31
      Abb33 = F1*F32
      Abb34 = F15*F29

      Abb14 = F19*F8
      Abb13 = F18*F6

      AbbQu4 = F19*F8
      AbbQl3 = F6*F18
      

      AbbSum52 = Abb23 + Abb25
      AbbSum14 = Abb23 + Abb29
      AbbSum45 = Abb33 - 1/2.D0*Abb29
      AbbSum181 = Abb25 - 1/2.D0*Abb34
      AbbSum88 = -Abb10 + Abb8
      AbbSum109 = -Abb24 + Abb8
      AbbSum74 = Abb8 - 2*Abb9
      AbbSum13 = Abb24 - Abb9
      AbbSum18 = Abb32 - Abb9
      AbbSum10 = Abb8 - Abb9
      AbbSum36 = Abb23 + Abb9
      AbbSum63 = Abb33 - AbbSum52
      AbbSum67 = -Abb34 + AbbSum14
      AbbSum32 = AbbSum52 - 1/2.D0*Abb34
      AbbSum184 = Abb23 - 1/3.D0*Abb29 + 1/3.D0*Abb34
      AbbSum175 = Abb23 - Abb29 + Abb34
      AbbSum162 = Abb32 - Abb33 - Abb8
      AbbSum107 = Abb8 - 1/3.D0*Abb10 - 2/3.D0*Abb24
      AbbSum98 = AbbSum109 + 1/2.D0*Abb25
      AbbSum95 = Abb25 + AbbSum109
      AbbSum4 = Abb10 + AbbSum74
      AbbSum78 = AbbSum13 - 1/2.D0*Abb10
      AbbSum17 = -Abb33 + AbbSum18
      AbbSum49 = -Abb29 + AbbSum10
      AbbSum164 = AbbSum10 + 1/3.D0*Abb29
      AbbSum154 = Abb29 + AbbSum10
      AbbSum173 = Abb33 + AbbSum10
      AbbSum102 = Abb8 - 1/3.D0*Abb10 - 2/3.D0*Abb9
      AbbSum105 = Abb8 - 1/2.D0*Abb24 - 1/2.D0*Abb9
      AbbSum56 = -Abb10 + AbbSum36
      AbbSum189 = -AbbSum181 + AbbSum45
      AbbSum42 = Abb32 - AbbSum45 - 3/2.D0*Abb9
      AbbSum167 = AbbSum10 + 1/2.D0*Abb23 + 1/2.D0*Abb25
      AbbSum172 = AbbSum10 + AbbSum52
      AbbSum111 = AbbSum109 + AbbSum18
      AbbSum116 = -Abb34 + AbbSum154
      AbbSum136 = AbbSum164 - 1/3.D0*Abb34
      AbbSum19 = Abb34 + AbbSum49
      AbbSum83 = -Abb24 + Abb32 + 1/2.D0*Abb10 - 1/2.D0*Abb9
      AbbSum11 = -AbbSum13 + AbbSum52
      AbbSum61 = -Abb24 + Abb32 - AbbSum63
      AbbSum126 = AbbSum109 + AbbSum181 + 1/2.D0*Abb10
      AbbSum171 = AbbSum164 + 2/3.D0*Abb23 + 2/3.D0*Abb25
      AbbSum149 = AbbSum10 + AbbSum32
      AbbSum142 = AbbSum167 - 1/4.D0*Abb34
      AbbSum97 = AbbSum105 + 1/2.D0*Abb23 + 1/2.D0*Abb25
      AbbSum161 = Abb25 + AbbSum105 + 1/2.D0*Abb23
      AbbSum93 = AbbSum107 + 1/3.D0*Abb23 + 2/3.D0*Abb25 - 1/3.D0*Abb29
      AbbSum180 = Abb23 - AbbSum181 + AbbSum78
      AbbSum146 = AbbSum171 - 2/3.D0*Abb34
      AbbSum29 = AbbSum32 - AbbSum78
      AbbSum178 = Abb10 - Abb23 - 2*Abb24 + 2*Abb25 + Abb29 - 
     -     2*Abb34 + Abb8
      AbbSum100 = AbbSum111 - AbbSum63
      AbbSum159 = Abb8 + 1/3.D0*Abb10 + 1/3.D0*Abb23 - 
     -     2/3.D0*Abb24 + 2/3.D0*Abb25 + 1/3.D0*Abb29 - 2/3.D0*Abb9
      AbbSum132 = AbbSum105 + AbbSum181 + 1/4.D0*Abb10 + 1/2.D0*Abb23
      AbbSum129 = AbbSum97 + 1/4.D0*Abb10 - 1/4.D0*Abb34
      AbbSum151 = -AbbSum149 + AbbSum45
      AbbSum123 = AbbSum159 - 2/3.D0*Abb34
      AbbSum193 = -AbbSum189 + AbbSum83
      AbbSum191 = Abb25 + Abb29 - Abb33 - Abb34 + AbbSum83 - 
     -     1/2.D0*Abb23
      AbbSum187 = -Abb23 + AbbSum193
      AbbSum72 = AbbSum32 - AbbSum45 + AbbSum83
      AbbSum182 = Abb10 - 2*Abb24 + 2*Abb25 + Abb32 - Abb34 + 
     -     Abb8 - AbbSum45 - 1/2.D0*Abb9
      AbbSum134 = AbbSum109 + AbbSum32 + AbbSum42 + 1/2.D0*Abb10

      AbbSumL13 = Abb10 - Abb34
      AbbSumL16 = -Abb29 + Abb9
      AbbSumL39 = AbbSumL13 + AbbSumL16
      AbbSumL12 = AbbSumL16 - 1/3.D0*Abb10 + 1/3.D0*Abb34
      AbbSumL21 = AbbSumL16 - 2/3.D0*Abb10 + 2/3.D0*Abb34
      AbbSumL3 = -AbbSumL13 + AbbSumL16

      AbbSumU27 = Abb10 - 2*Abb34
      AbbSumU8 = Abb10 - Abb34
      AbbSumU34 = Abb10 - 1/2.D0*Abb34
      AbbSumU12 = Abb10 + 1/2.D0*Abb34
      AbbSumU22 = Abb10 + 2*Abb34


      end


*******************************************************************************
*******************************************************************************

*** Subroutine that calculates the spinor products for qbar + q -> qbar + H + q
*** (qk = 2)

      subroutine hel_qbarqCC

      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel

** Vector products  
      double complex F4, F17, F22, F6, F7, F20, F23, F19, F2, F21
      double complex F10, F34, F13, F15, F11, F5, F12, F18, F1, F16
      double complex F8, F31, F24, F25, F29, F3, F30, F14, F32, F33
      double complex F9, F26, F27, F28    

** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS


** Contains the declarations of the abbreviations
#include "abbr_CC.h"


      F2 = Hel(1)*vSeS(DottedSpinor(1,-1,7),
     -     DottedSpinor(5,1,7))
      F1 = Hel(3)*vSeS(Spinor(3,-1,7),Spinor(2,1,7))

      F4 = Hel(1)*vSxS(DottedSpinor(1,-1,7),Spinor(2,1,6))
      F22 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -     Spinor(2,1,7)))
      F20 = -(Hel(3)*vSxS(DottedSpinor(5,1,7),
     -     Spinor(3,-1,6)))
      F19 = Hel(3)*vSxS(DottedSpinor(5,1,6),Spinor(3,-1,7))
      F10 = Hel(1)*Hel(3)*
     -     vSxS(DottedSpinor(1,-1,7),vVxS(k(2),
     -     Spinor(3,-1,7)))
      F34 = Hel(1)*vSxS(DottedSpinor(1,-1,7),vVxS(k(4),
     -     Spinor(2,1,7)))
      F13 = Hel(1)*Hel(3)*
     -     vSxS(DottedSpinor(1,-1,7),vVxS(k(4),
     -     Spinor(3,-1,7)))
      F11 = vSxS(DottedSpinor(5,1,7),vVxS(k(1),
     -     Spinor(2,1,7)))
      F5 = Hel(3)*vSxS(DottedSpinor(5,1,7),vVxS(k(1),
     -     Spinor(3,-1,7)))
      F12 = vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(2,1,7)))
      F18 = Hel(3)*vSxS(DottedSpinor(5,1,7),vVxS(k(3),
     -     Spinor(3,-1,7)))
      F16 = -(Hel(3)*vSeS(Spinor(3,-1,6),Spinor(2,1,6)))
      F31 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVeS(k(2),DottedSpinor(5,1,6)))
      F24 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVeS(k(3),DottedSpinor(5,1,6)))
      F25 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVeS(k(4),DottedSpinor(5,1,6)))
      F29 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -     vBeS(k(1),DottedSpinor(5,1,7))))
      F3 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -     vBeS(k(2),DottedSpinor(5,1,7))))
      F30 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -     vBeS(k(3),DottedSpinor(5,1,7))))
      F14 = -(Hel(1)*vSxS(DottedSpinor(1,-1,6),
     -     vBeS(k(4),DottedSpinor(5,1,7))))
      F32 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(1),DottedSpinor(5,1,7))))
      F33 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(3),DottedSpinor(5,1,7))))
      F9 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(4),DottedSpinor(5,1,7))))
      F26 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(3),vBeS(k(1),DottedSpinor(5,1,7))))
      F27 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),vBeS(k(1),DottedSpinor(5,1,7))))
      F28 = Hel(1)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),vBeS(k(3),DottedSpinor(5,1,7))))


      Abb1 = F1*F2
      Abb37 = F1*F32
      Abb35 = F11*F13
      Abb34 = F1*F9
      Abb31 = F1*F28
      Abb29 = F1*F26
      Abb27 = F1*F33
      Abb10 = F1*F27
      Abb9 = F10*F11
      Abb5 = F10*F12
      Abb2 = F12*F13

      Abb21 = F20*F22
      Abb8 = F19*F4

      AbbQu3 = F22*F20
      AbbQl3 = F4*F19


      AbbSum79 = Abb29 + Abb10
      AbbSum30 = Abb37 - 1/2.D0*Abb9
      AbbSum43 = Abb29 + Abb9
      AbbSum42 = Abb31 - Abb5
      AbbSum46 = Abb27 - Abb5
      AbbSum20 = Abb29 + Abb5
      AbbSum127 = Abb2 + Abb34
      AbbSum149 = Abb31 + Abb34
      AbbSum165 = Abb34 + 2/3.D0*Abb5
      AbbSum9 = Abb5 + Abb34
      AbbSum48 = 2*Abb5 + Abb34
      AbbSum83 = Abb37 - AbbSum79
      AbbSum190 = -Abb29 - Abb31 + Abb27
      AbbSum178 = Abb2 + Abb29 - Abb5
      AbbSum53 = AbbSum42 - 1/2.D0*Abb2
      AbbSum184 = Abb29 + AbbSum42
      AbbSum187 = Abb29 + 1/3.D0*Abb2 - 1/3.D0*Abb5
      AbbSum69 = -Abb2 + AbbSum20
      AbbSum45 = Abb37 - AbbSum46
      AbbSum59 = -Abb35 + AbbSum43
      AbbSum16 = AbbSum79 - 1/2.D0*Abb35
      AbbSum171 = Abb34 - 1/3.D0*Abb2 + 2/3.D0*Abb31
      AbbSum147 = Abb34 + 1/3.D0*Abb2 + 2/3.D0*Abb31
      AbbSum173 = AbbSum149 - 1/2.D0*Abb2
      AbbSum135 = -Abb10 + AbbSum149
      AbbSum138 = AbbSum149 - 1/2.D0*Abb10
      AbbSum175 = Abb27 + Abb34 - 1/2.D0*Abb5
      AbbSum145 = Abb34 + 1/2.D0*Abb31 + 1/2.D0*Abb5
      AbbSum141 = AbbSum165 + 1/3.D0*Abb2
      AbbSum76 = Abb9 + AbbSum9
      AbbSum33 = -Abb2 + AbbSum48
      AbbSum194 = -Abb10 + Abb31 + Abb37 - Abb27
      AbbSum183 = -Abb10 + AbbSum184
      AbbSum65 = -Abb31 + Abb27 + 1/2.D0*Abb2 - 1/2.D0*Abb5
      AbbSum40 = -AbbSum42 + AbbSum79
      AbbSum27 = -Abb27 + AbbSum30 + 3/2.D0*Abb5
      AbbSum151 = -AbbSum149 + AbbSum46
      AbbSum168 = AbbSum145 - 1/4.D0*Abb2
      AbbSum1 = -Abb35 + AbbSum76
      AbbSum109 = AbbSum9 - 1/3.D0*Abb9 + 1/3.D0*Abb35
      AbbSum89 = -Abb9 + Abb35 + AbbSum9
      AbbSum74 = Abb31 - Abb27 + AbbSum83
      AbbSum189 = -Abb10 + Abb37 - AbbSum190
      AbbSum163 = -AbbSum175 + AbbSum30
      AbbSum185 = -2*Abb10 + 2*Abb31 + Abb37 - Abb27 + Abb34
      AbbSum136 = AbbSum145 - 1/2.D0*Abb29 - 1/2.D0*Abb10
      AbbSum162 = AbbSum173 - 1/2.D0*Abb10 + 1/4.D0*Abb35
      AbbSum115 = AbbSum9 - 1/2.D0*Abb29 - 1/2.D0*Abb10 + 1/4.D0*Abb35
      AbbSum99 = -Abb10 + AbbSum173 + 1/2.D0*Abb35
      AbbSum122 = -AbbSum16 + AbbSum9
      AbbSum13 = AbbSum16 - AbbSum53
      AbbSum159 = AbbSum171 - 1/3.D0*Abb29 - 2/3.D0*Abb10 + 1/3.D0*Abb9
      AbbSum132 = AbbSum147 - 1/3.D0*Abb29 - 2/3.D0*Abb10 + 1/3.D0*Abb9
      AbbSum119 = AbbSum9 - 2/3.D0*Abb29 - 2/3.D0*Abb10 - 
     -     1/3.D0*Abb9 + 2/3.D0*Abb35
      AbbSum192 = AbbSum194 + 1/2.D0*Abb2 + 1/2.D0*Abb29 - 1/2.D0*Abb5
      AbbSum180 = -2*Abb10 + 2*Abb31 - 2*Abb5 + AbbSum127 + AbbSum43
      AbbSum139 = -AbbSum151 + AbbSum83
      AbbSum101 = AbbSum168 - 1/2.D0*Abb29 - 1/2.D0*Abb10 + 1/4.D0*Abb35
      AbbSum105 = -Abb10 + AbbSum168 - 1/2.D0*Abb29 + 1/2.D0*Abb35
      AbbSum124 = AbbSum122 + AbbSum30
      AbbSum96 = AbbSum165 - 1/3.D0*Abb2 - 1/3.D0*Abb29 - 
     -     2/3.D0*Abb10 + 2/3.D0*Abb31 - 1/3.D0*Abb9 + 2/3.D0*Abb35
      AbbSum63 = -AbbSum16 + AbbSum30 - AbbSum65
      AbbSum106 = -AbbSum16 + AbbSum173 + AbbSum27

      AbbSumL13 = -Abb2 + Abb35
      AbbSumL16 = Abb9 - Abb5
      AbbSumL3 = -AbbSumL13 + AbbSumL16
      AbbSumL21 = AbbSumL16 + 2/3.D0*Abb2 - 2/3.D0*Abb35
      AbbSumL12 = AbbSumL16 + 1/3.D0*Abb2 - 1/3.D0*Abb35
      AbbSumL39 = AbbSumL13 + AbbSumL16

      AbbSumU26 = Abb35 - 2*Abb2
      AbbSumU7 = Abb35 - Abb2
      AbbSumU34 = Abb35 - 1/2.D0*Abb2
      AbbSumU11 = Abb35 + 1/2.D0*Abb2
      AbbSumU22 = Abb35 + 2*Abb2

      
      end


******************************************************************************
******************************************************************************

*** Subroutine calculating spinor products for q + qbar -> q + H + qbar
*** (qk = 3)

      subroutine hel_qqbarCC

      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel

** Vector products
      double complex F4, F17, F22, F6, F7, F20, F9, F19, F2, F21
      double complex F10, F34, F26, F15, F12, F5, F13, F18, F1, F16
      double complex F27, F8, F28, F29, F3, F30, F31, F14, F32, F23
      double complex F24, F33, F25, F11
      
** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS


** Contains the declarations of the abbreviations
#include "abbr_CC.h"


      F2 = Hel(2)*vSeS(DottedSpinor(2,-1,7),
     -     DottedSpinor(3,1,7))
      F1 = Hel(5)*vSeS(Spinor(5,-1,7),Spinor(1,1,7))

      F4 = Hel(2)*vSxS(DottedSpinor(2,-1,7),Spinor(1,1,6))
      F22 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -     Spinor(1,1,7)))
      F20 = -(Hel(5)*vSxS(DottedSpinor(3,1,7),
     -     Spinor(5,-1,6)))
      F19 = Hel(5)*vSxS(DottedSpinor(3,1,6),Spinor(5,-1,7))
      F10 = Hel(2)*Hel(5)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(1),Spinor(5,-1,7)))
      F34 = Hel(2)*vSxS(DottedSpinor(2,-1,7),vVxS(k(3),
     -     Spinor(1,1,7)))
      F26 = Hel(2)*Hel(5)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(3),Spinor(5,-1,7)))
      F15 = -(Hel(2)*vSeS(DottedSpinor(2,-1,6),
     -     DottedSpinor(3,1,6)))
      F12 = vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     Spinor(1,1,7)))
      F5 = Hel(5)*vSxS(DottedSpinor(3,1,7),vVxS(k(2),
     -     Spinor(5,-1,7)))
      F13 = vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(1,1,7)))
      F18 = Hel(5)*vSxS(DottedSpinor(3,1,7),vVxS(k(4),
     -     Spinor(5,-1,7)))
      F16 = -(Hel(5)*vSeS(Spinor(5,-1,6),Spinor(1,1,6)))
      F27 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVeS(k(1),DottedSpinor(3,1,6)))
      F28 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVeS(k(3),DottedSpinor(3,1,6)))
      F29 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVeS(k(4),DottedSpinor(3,1,6)))
      F3 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -     vBeS(k(1),DottedSpinor(3,1,7))))
      F31 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -     vBeS(k(3),DottedSpinor(3,1,7))))
      F14 = -(Hel(2)*vSxS(DottedSpinor(2,-1,6),
     -     vBeS(k(4),DottedSpinor(3,1,7))))
      F32 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(1),vBeS(k(2),DottedSpinor(3,1,7))))
      F23 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(1),vBeS(k(3),DottedSpinor(3,1,7))))
      F24 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(1),vBeS(k(4),DottedSpinor(3,1,7))))
      F33 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(3),vBeS(k(2),DottedSpinor(3,1,7))))
      F25 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(3),vBeS(k(4),DottedSpinor(3,1,7))))
      F11 = Hel(2)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(4),vBeS(k(2),DottedSpinor(3,1,7))))


      Abb1 = F1*F2

      Abb27 = F20*F22
      Abb14 = F19*F4

      Abb34 = F13*F26
      Abb33 = F1*F32
      Abb26 = F1*F24
      Abb24 = F1*F33
      Abb16 = F1*F25
      Abb11 = F1*F23
      Abb10 = F10*F13
      Abb8 = F10*F12
      Abb7 = F1*F11
      Abb5 = F12*F26

      AbbQu4 = F4*F19
      AbbQl3 = F22*F20


      Opt1 = -Abb26 - 2*Abb7 + Abb8 + Abb10
      Opt2 = -2*Abb26 + 2*Abb16 - 3*Abb7 + 2*Abb8 + 2*Abb10
      AbbSum44 = Abb11 + Abb16
      AbbSum22 = Abb11 - Abb5
      AbbSum52 = Abb33 - Abb24
      AbbSum210 = -Abb16 + Abb24
      AbbSum80 = Abb7 - 2*Abb8
      AbbSum47 = Abb11 - Abb8
      AbbSum11 = Abb7 - Abb8
      AbbSum30 = Abb33 - 1/2.D0*Abb8
      AbbSum90 = Abb7 - Abb10
      AbbSum60 = AbbSum44 - 1/2.D0*Abb34
      AbbSum76 = Abb34 + AbbSum22
      AbbSum206 = Abb11 - 1/3.D0*Abb34 + 1/3.D0*Abb5
      AbbSum191 = Abb11 - Abb34 + Abb5
      AbbSum57 = -Abb5 + AbbSum11
      AbbSum176 = AbbSum11 + 1/3.D0*Abb5
      AbbSum164 = Abb5 + AbbSum11
      AbbSum120 = Abb26 - 2*Abb10 + 3/2.D0*Abb7
      AbbSum122 = Abb26 + Abb7 - 3/2.D0*Abb10
      AbbSum18 = Abb26 + Abb8 - 3/2.D0*Abb10
      AbbSum84 = Abb26 + Abb8 - Abb10
      AbbSum113 = Abb7 - 2/3.D0*Abb8 - 1/3.D0*Abb10
      AbbSum2 = Abb10 + AbbSum80
      AbbSum64 = Abb10 + AbbSum47
      AbbSum34 = Abb34 + AbbSum57
      AbbSum149 = AbbSum176 - 1/3.D0*Abb34
      AbbSum129 = -Abb34 + AbbSum164
      AbbSum29 = AbbSum52 + 1/2.D0*Abb5 - 1/2.D0*Abb8
      AbbSum159 = -AbbSum11 + AbbSum44
      AbbSum188 = Abb24 - AbbSum11 - 1/2.D0*Abb5
      AbbSum153 = -2*Abb7 + 2*Abb8 + AbbSum44
      AbbSum117 = Abb26 + 2*Abb7 - Abb8 - 3/2.D0*Abb10
      AbbSum140 = Abb26 - Abb16 + AbbSum90
      AbbSum87 = Abb33 - AbbSum84
      AbbSum73 = -Abb24 + AbbSum60 + 1/2.D0*Abb5
      AbbSum173 = Abb7 + AbbSum52 + 1/2.D0*Abb5 - 3/2.D0*Abb8
      AbbSum186 = -AbbSum11 + AbbSum60
      AbbSum161 = -Abb24 + AbbSum159
      AbbSum179 = -2*Abb7 + 2*Abb8 + AbbSum60
      AbbSum106 = AbbSum122 - 1/2.D0*Abb16 + 1/4.D0*Abb34
      AbbSum99 = -Abb16 + AbbSum122 + 1/2.D0*Abb34
      AbbSum199 = Abb11 - Abb16 + AbbSum84
      AbbSum41 = AbbSum44 - AbbSum84
      AbbSum124 = -AbbSum122 + AbbSum30
      AbbSum183 = AbbSum60 - 1/2.D0*Abb5 - 3/2.D0*Abb7 + 3/2.D0*Abb8
      AbbSum156 = AbbSum44 + 1/2.D0*Abb34 - 1/2.D0*Abb5 - 
     -     3/2.D0*Abb7 + 3/2.D0*Abb8
      AbbSum215 = Abb16 + AbbSum52 - AbbSum84
      AbbSum141 = AbbSum44 + Opt1
      AbbSum144 = Abb11 + 2*Abb16 + Opt1
      AbbSum15 = -AbbSum18 + AbbSum60
      AbbSum195 = 2*Abb26 - 2*Abb16 + Abb7 - 2*Abb10 + AbbSum191
      AbbSum146 = -Abb26 + AbbSum44 + AbbSum52 - AbbSum90
      AbbSum208 = -AbbSum199 + AbbSum52
      AbbSum79 = AbbSum41 + AbbSum52
      AbbSum172 = Abb11 + 2*Abb16 - Abb34 - AbbSum117
      AbbSum101 = -AbbSum117 + AbbSum60
      AbbSum202 = -2*Abb26 + 2*Abb16 - Abb7 - Abb8 + 2*Abb10 + AbbSum52
      AbbSum94 = -2*Abb26 + 2*Abb16 - 3*Abb7 + 4*Abb10 + AbbSum191
      AbbSum168 = -Abb34 + AbbSum22 + Opt2
      AbbSum136 = AbbSum76 + Opt2
      AbbSum212 = AbbSum215 - 1/2.D0*Abb11 + 1/2.D0*Abb34 - 1/2.D0*Abb5
      AbbSum70 = -Abb26 + AbbSum52 + AbbSum60 + 1/2.D0*Abb5 - 
     -     3/2.D0*Abb8 + 3/2.D0*Abb10
      AbbSum108 = -AbbSum122 + AbbSum29 + AbbSum60

      AbbSumL16 = Abb8 - Abb5
      AbbSumL13 = Abb10 - Abb34
      AbbSumL39 = AbbSumL13 + AbbSumL16
      AbbSumL12 = AbbSumL16 - 1/3.D0*Abb10 + 1/3.D0*Abb34
      AbbSumL21 = AbbSumL16 - 2/3.D0*Abb10 + 2/3.D0*Abb34
      AbbSumL1 = -AbbSumL13 + AbbSumL16

      AbbSumU27 = Abb10 - 2*Abb34
      AbbSumU8 = Abb10 - Abb34
      AbbSumU34 = Abb10 - 1/2.D0*Abb34
      AbbSumU12 = Abb10 + 1/2.D0*Abb34
      AbbSumU22 = Abb10 + 2*Abb34


      end


******************************************************************************
******************************************************************************

*** Subroutine to calculate spinor products for qbar + qbar -> qbar + H + qbar
*** (qk = 4)

      subroutine hel_qbarqbarCC

      implicit none


** Helicities of the particles
      integer Hel(5)
      common /FChelicities/ Hel

** Vector products
      double complex F16, F19, F15, F22, F17, F24, F21, F8, F2, F9
      double complex F18, F10, F31, F5, F20, F12, F25, F13, F1, F23
      double complex F3, F32, F6, F7, F29, F4, F30, F14, F33, F34
      double complex F11, F26, F27, F28

** The spinors' and dotted spinors' id's
      integer DottedSpinor, Spinor, k
      external DottedSpinor
      external Spinor
      external k

** The vector products
      double complex vSeS, vSxS
      integer vVxS, vVeS, vBeS
      external vSeS
      external vSxS
      external vVxS
      external vVeS
      external vBeS


** Contains the declarations of the abbreviations
#include "abbr_CC.h"


      F2 = Hel(1)*Hel(2)*
     -     vSeS(DottedSpinor(1,-1,7),DottedSpinor(2,-1,7))
      F1 = Hel(3)*Hel(5)*vSeS(Spinor(3,-1,7),Spinor(5,-1,7))

      F19 = -(Hel(1)*Hel(5)*vSxS(DottedSpinor(1,-1,7),
     -     Spinor(5,-1,6)))
      F22 = -(Hel(1)*Hel(5)*vSxS(DottedSpinor(1,-1,6),
     -     Spinor(5,-1,7)))
      F17 = -(Hel(2)*Hel(3)*vSxS(DottedSpinor(2,-1,7),
     -     Spinor(3,-1,6)))
      F21 = -(Hel(2)*Hel(3)*vSxS(DottedSpinor(2,-1,6),
     -     Spinor(3,-1,7)))
      F9 = Hel(1)*Hel(3)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),Spinor(3,-1,7)))
      F18 = Hel(1)*Hel(5)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),Spinor(5,-1,7)))
      F10 = Hel(1)*Hel(3)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),Spinor(3,-1,7)))
      F31 = Hel(1)*Hel(5)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),Spinor(5,-1,7)))
      F5 = Hel(1)*Hel(2)*
     -     vSeS(DottedSpinor(1,-1,6),DottedSpinor(2,-1,6))
      F20 = Hel(2)*Hel(3)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(1),Spinor(3,-1,7)))
      F12 = Hel(2)*Hel(5)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(1),Spinor(5,-1,7)))
      F25 = Hel(2)*Hel(3)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(3),Spinor(3,-1,7)))
      F13 = Hel(2)*Hel(5)*vSxS(DottedSpinor(2,-1,7),
     -     vVxS(k(3),Spinor(5,-1,7)))
      F23 = Hel(3)*Hel(5)*vSeS(Spinor(3,-1,6),Spinor(5,-1,6))
      F6 = -(Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVeS(k(3),DottedSpinor(2,-1,6))))
      F7 = -(Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVeS(k(4),DottedSpinor(2,-1,6))))
      F30 = -(Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,6),
     -     vBeS(k(3),DottedSpinor(2,-1,7))))
      F14 = -(Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,6),
     -     vBeS(k(4),DottedSpinor(2,-1,7))))
      F33 = Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(1),DottedSpinor(2,-1,7))))
      F34 = Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(3),DottedSpinor(2,-1,7))))
      F11 = Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(2),vBeS(k(4),DottedSpinor(2,-1,7))))
      F26 = Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(3),vBeS(k(1),DottedSpinor(2,-1,7))))
      F27 = Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),vBeS(k(1),DottedSpinor(2,-1,7))))
      F28 = Hel(1)*Hel(2)*vSxS(DottedSpinor(1,-1,7),
     -     vVxS(k(4),vBeS(k(3),DottedSpinor(2,-1,7))))


      Abb1 = F1*F2

      Abb36 = F12*F9
      Abb35 = F10*F12
      Abb34 = F1*F33
      Abb28 = F1*F27
      Abb22 = F1*F26
      Abb15 = F13*F9
      Abb10 = F1*F34
      Abb9 = F1*F11
      Abb8 = F1*F28
      Abb2 = F10*F13

      Abb27 = F19*F21
      Abb21 = F17*F22

      AbbQu3 = F22*F17
      AbbQl3 = F19*F21


      Opt1 = -(1/3.D0*Abb35) + 1/3.D0*Abb36 + 1/3.D0*Abb22 + 2/3.D0*Abb28
      AbbSum45 = -Abb36 + Abb22
      AbbSum19 = -Abb15 + Abb22
      AbbSum15 = Abb22 + Abb28
      AbbSum24 = Abb34 - Abb10
      AbbSum52 = Abb10 - 1/2.D0*Abb15
      AbbSum175 = Abb9 + 2/3.D0*Abb15
      AbbSum9 = Abb15 + Abb9
      AbbSum57 = 2*Abb15 + Abb9
      AbbSum125 = Abb2 + Abb9
      AbbSum64 = Abb35 + AbbSum45
      AbbSum187 = Abb15 - Abb2 + Abb22
      AbbSum204 = Abb22 + 1/3.D0*Abb15 - 1/3.D0*Abb2
      AbbSum71 = Abb2 + AbbSum19
      AbbSum81 = AbbSum15 - 1/2.D0*Abb35
      AbbSum41 = Abb15 + Abb8 - 3/2.D0*Abb2
      AbbSum61 = Abb15 - Abb2 + Abb8
      AbbSum184 = Abb10 - AbbSum9
      AbbSum77 = Abb36 + AbbSum9
      AbbSum27 = -Abb2 + AbbSum57
      AbbSum148 = AbbSum175 + 1/3.D0*Abb2
      AbbSum182 = -Abb8 + AbbSum125
      AbbSum157 = -Abb8 + Abb9 + 3/2.D0*Abb2
      AbbSum180 = Abb9 + 2/3.D0*Abb2 - 2/3.D0*Abb8
      AbbSum155 = Abb9 + 4/3.D0*Abb2 - 2/3.D0*Abb8
      AbbSum198 = Abb22 + AbbSum41
      AbbSum51 = AbbSum24 - 1/2.D0*Abb36 + 1/2.D0*Abb15
      AbbSum69 = Abb10 - AbbSum61
      AbbSum89 = Abb35 - Abb36 + AbbSum9
      AbbSum110 = AbbSum9 + 1/3.D0*Abb35 - 1/3.D0*Abb36
      AbbSum2 = -Abb35 + AbbSum77
      AbbSum114 = AbbSum9 + 1/2.D0*Abb22 + 1/2.D0*Abb28
      AbbSum120 = AbbSum15 + AbbSum9
      AbbSum171 = AbbSum182 + 1/2.D0*Abb28
      AbbSum100 = Abb28 + AbbSum182
      AbbSum178 = Abb9 + 1/2.D0*Abb15 + 1/2.D0*Abb2 - 1/2.D0*Abb8
      AbbSum152 = Abb9 + 1/2.D0*Abb15 + 3/4.D0*Abb2 - 1/2.D0*Abb8
      AbbSum173 = AbbSum24 + AbbSum9
      AbbSum12 = Abb28 - Abb8 + AbbSum71
      AbbSum87 = Abb34 - AbbSum81 - 1/2.D0*Abb36
      AbbSum209 = -Abb22 - Abb8 + Abb10 - 3/2.D0*Abb15 + 3/2.D0*Abb2
      AbbSum122 = Abb34 - AbbSum120
      AbbSum142 = AbbSum157 - 1/4.D0*Abb35 + 1/2.D0*Abb28
      AbbSum135 = Abb28 + AbbSum157 - 1/2.D0*Abb35
      AbbSum159 = AbbSum157 + AbbSum52
      AbbSum38 = Abb28 - Abb8 + AbbSum19 - 1/2.D0*Abb35 + 3/2.D0*Abb2
      AbbSum195 = -Abb28 + AbbSum198 + 1/2.D0*Abb35
      AbbSum117 = AbbSum110 + 2/3.D0*Abb22 + 2/3.D0*Abb28
      AbbSum101 = AbbSum178 + 1/2.D0*Abb22 + 1/2.D0*Abb28
      AbbSum105 = Abb28 + AbbSum178 + 1/2.D0*Abb22
      AbbSum68 = -AbbSum12 + AbbSum24
      AbbSum107 = -AbbSum15 - AbbSum182 + AbbSum24
      AbbSum167 = AbbSum180 + Opt1
      AbbSum131 = AbbSum155 + Opt1
      AbbSum137 = AbbSum152 - 1/4.D0*Abb35 + 1/2.D0*Abb22 + 1/2.D0*Abb28
      AbbSum213 = -Abb28 + Abb8 + AbbSum24 + 1/2.D0*Abb35 - 
     -     1/2.D0*Abb36 + 3/2.D0*Abb15 - 3/2.D0*Abb2
      AbbSum191 = -2*Abb15 + 4*Abb2 + 2*Abb28 - 2*Abb8 + 
     -     Abb9 - AbbSum64
      AbbSum96 = AbbSum175 + 1/3.D0*Abb35 - 1/3.D0*Abb36 + 
     -     2/3.D0*Abb2 + 1/3.D0*Abb22 + 2/3.D0*Abb28 - 2/3.D0*Abb8
      AbbSum76 = Abb8 + AbbSum24 - AbbSum81 - 1/2.D0*Abb36 + 
     -     3/2.D0*Abb15 - 3/2.D0*Abb2
      AbbSum210 = 2*Abb15 - 2*Abb2 - Abb28 + Abb8 + 
     -     AbbSum24 + 1/2.D0*Abb35 - 1/2.D0*Abb36 + 1/2.D0*Abb22
      AbbSum207 = Abb22 + AbbSum213
      AbbSum199 = Abb35 - 3*Abb2 - 2*Abb28 + 2*Abb8 - Abb9 + 
     -     AbbSum24 - 1/2.D0*Abb36 + 3/2.D0*Abb15
      AbbSum144 = -AbbSum157 + AbbSum51 - AbbSum81

      AbbSumL16 = -Abb15 + Abb36
      AbbSumL13 = -Abb2 + Abb35
      AbbSumL1 = -AbbSumL13 + AbbSumL16
      AbbSumL21 = AbbSumL16 + 2/3.D0*Abb2 - 2/3.D0*Abb35
      AbbSumL12 = AbbSumL16 + 1/3.D0*Abb2 - 1/3.D0*Abb35
      AbbSumL39 = AbbSumL13 + AbbSumL16

      AbbSumU26 = Abb35 - 2*Abb2
      AbbSumU7 = Abb35 - Abb2
      AbbSumU34 = Abb35 - 1/2.D0*Abb2
      AbbSumU11 = Abb35 + 1/2.D0*Abb2
      AbbSumU22 = Abb35 + 2*Abb2


      end


*******************************************************************************
