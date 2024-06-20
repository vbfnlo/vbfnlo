********************************************************************************

********************************************************************
*** ew_virtual.F                                                 ***
*** sophy@particle.uni-karlsruhe.de                              ***
***                                                              ***
*** This set of subroutines pulls together all of the virtual    ***
*** corrections into a nice, big file.                           ***
********************************************************************

********************************************************************************

*** This subroutine calls and combines all contributions to the HVV form factors

      subroutine form_factors(q1, q2, ahvvL, ewNLO)

      implicit none


** formfactors.  Labelling: (i,2,2) are the z formfactors, (i,3,4) and (i,4,3)
** are the w formfactors.
**    i     =     1      a1
**          =     2      a2
**          =     3      a3
      double complex ahvvL(3,4,4)


** These are the contributions to the formfactors from the tree level (the 
** first entry in each array is for the light higgs, the second is for the 
** heavy higgs) for W and Z bosons respectively.  
      double complex treeW(3), treeZ(3)

** These are the contributions to the formfactors from the vertex, for w fusion,
** zfusion, photon fusion and photon-z fusion
      double complex a1vtxW, a2vtxW, a3vtxW
      double complex a1vtxZ, a2vtxZ, a3vtxZ
      double complex a1vtxA, a2vtxA, a3vtxA
      double complex a1vtxAZ, a2vtxAZ, a3vtxAZ
      double complex a1vtxAZ2, a2vtxAZ2, a3vtxAZ2
      double complex a1zh, a1gamh

** These are the contributions to the formfactors from the vector self energy, 
** for w fusion, z fusion and photon-z fusion
      double complex a1slfW, a1slfZ, a1slfAZ, a1slfAZ2

** These are the contributions to the formfactors from the higgs propagator 
** factors, for w fusion and zfusion
      double complex a1propW, a1propZ


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions
      double complex ctFM
      double complex a1vtxWfm, a2vtxWfm, a3vtxWfm
      double complex a1vtxZfm, a2vtxZfm, a3vtxZfm
      double complex a1vtxAHfm, a2vtxAHfm, a3vtxAHfm
      double complex a1vtxAZHfm, a2vtxAZHfm, a3vtxAZHfm
      double complex a1vtxAZH2fm, a2vtxAZH2fm, a3vtxAZH2fm
      double complex a1selfWfm, a1selfZfm, a1selfAZfm, a1selfAZ2fm
      common /fmloop/ a1vtxWfm, a2vtxWfm, a3vtxWfm
      common /fmloop/ a1vtxZfm, a2vtxZfm, a3vtxZfm
      common /fmloop/ a1vtxAHfm, a2vtxAHfm, a3vtxAHfm
      common /fmloop/ a1vtxAZHfm, a2vtxAZHfm, a3vtxAZHfm
      common /fmloop/ a1vtxAZH2fm, a2vtxAZH2fm, a3vtxAZH2fm
      common /fmloop/ a1selfWfm, a1selfZfm, a1selfAZfm, a1selfAZ2fm
** formfactors for FERMIONS.  Labelling: (i,2,2) are the z formfactors, 
** (i,3,4) and (i,4,3) are the w formfactors.
**    i     =     1      a1
**          =     2      a2
**          =     3      a3
      double complex ahvvFM(3,4,4)  
      common /fmFormFac/ ahvvFM


** The momenta: q1 and q2 are the momenta of the vector bosons, q3 the momentum
** of the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2.  t is a do-loop variable.
      double precision q1(0:4), q2(0:4), q3(0:4), q12, q22, q32, q1q2
      integer t


** Loop parameter
      integer born

** Parameter controlling whether or not the virtual ew corrections are
** included
      integer ewNLO

** Debugging flag
      integer debug


** These are parameter declarations. 
#include "mssm.inc"
#include "global.inc"


      debug = 0


** Clearing the values of the tree level piece, the form factors, the VV-higgs
** vertex and the higgs propagator contribution
      do t = 1, 3
         treeW(t) = (0D0,0D0)
         treeZ(t) = (0D0,0D0)
      end do
      
      a1vtxW = (0D0,0D0)
      a1vtxZ = (0D0,0D0)
      a1vtxA = (0D0,0D0)
      a1vtxAZ = (0D0,0D0)
      a1vtxAZ2 = (0D0,0D0)
      a2vtxW = (0D0,0D0)
      a2vtxZ = (0D0,0D0)
      a2vtxA = (0D0,0D0)
      a2vtxAZ = (0D0,0D0)
      a2vtxAZ2 = (0D0,0D0)
      a3vtxW = (0D0,0D0)
      a3vtxZ = (0D0,0D0)
      a3vtxA = (0D0,0D0)
      a3vtxAZ = (0D0,0D0)
      a3vtxAZ2 = (0D0,0D0)

      a1slfW = (0D0,0D0)
      a1slfZ = (0D0,0D0)
      a1slfAZ = (0D0,0D0)
      a1slfAZ2 = (0D0,0D0)

      a1zh = (0D0,0D0)
      a1gamh = (0D0,0D0)
      
      a1propW = (0D0,0D0)
      a1propZ = (0D0,0D0)


** The momenta (q1,q2: w-bosons, q3:h-boson)
      do t = 0, 3
         q3(t) = q1(t) + q2(t)
      end do

      q12 = q1(0)**2
      q22 = q2(0)**2
      q32 = q3(0)**2
      q1q2 = q1(0)*q2(0)
      
      do t = 1, 3
         q12 = q12-(q1(t)**2)
         q22 = q22-(q2(t)**2)
         q32 = q32-(q3(t)**2)
         q1q2 = q1q2-(q1(t)*q2(t))
      end do
      
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'q12, q22, q32, q1q2 =', q12, q22, q32, q1q2
      end if


** Setting to tree level Higgs masses
      if ((model .eq. 2) .and. (MHloops .eq. 0)) then  ! working in the MSSM
         do t = 1, 4
            MHiggs(t) = MHiggs_tree(t)
         end do

         Mh0 = MHiggs(1)
         MHH = MHiggs(2)
         MA0 = MHiggs(3)
         MHp = MHiggs(4)
         
         Mh02 = Mh0**2
         MHH2 = MHH**2
         MA02 = MA0**2
         MHp2 = MHp**2
      end if
** End of setting tree level Higgs masses     

      
** For the standard model:
      if (model .eq. 1) then
         
         if (HiggsType .eq. 0) then
            call treeSM(treeW,treeZ,debug,1)
         else
            call treeSM(treeW,treeZ,debug,HiggsType)
         end if

         if (ewNLO .ne. 0) then
            call ww_vtxSM(a1vtxW,a2vtxW,a3vtxW,q12,q22,q32,q1q2,debug)
            call zz_vtxSM(a1vtxZ,a2vtxZ,a3vtxZ,q12,q22,q32,q1q2,debug)
            call aa_vtxSM(a1vtxA,a2vtxA,a3vtxA,q12,q22,q32,q1q2,debug)
            call az_vtxSM(a1vtxAZ,a2vtxAZ,a3vtxAZ,q12,q22,q32,q1q2,
     &           debug,1)
            call az_vtxSM(a1vtxAZ2,a2vtxAZ2,a3vtxAZ2,q22,q12,q32,q1q2,
     &           debug,2)

            call ww_selfSM(a1slfW,q12,q22,q32,q1q2,debug)
            call zz_selfSM(a1slfZ,q12,q22,q32,q1q2,debug)
            call az_selfSM(a1slfAZ,q12,q22,q32,q1q2,debug,1)
            call az_selfSM(a1slfAZ2,q22,q12,q32,q1q2,debug,2)
         end if
            
** For susy:
      else if (model .eq. 2) then

         call treeSUSY(treeW,treeZ,debug,HiggsType) 

         if (propLevel .ne. 0) then
            call hgs_prop(a1propW,a1propZ,treeW,treeZ,Hgs,debug,
     -           HiggsType)
         end if

         if (ewNLO .ne. 0) then
            call ww_vtxSUSY(a1vtxW,a2vtxW,a3vtxW,q12,q22,q32,q1q2,debug)
            call zz_vtxSUSY(a1vtxZ,a2vtxZ,a3vtxZ,q12,q22,q32,q1q2,debug)      
            call aa_vtxSUSY(a1vtxA,a2vtxA,a3vtxA,q12,q22,q32,q1q2,debug)
            call az_vtxSUSY(a1vtxAZ,a2vtxAZ,a3vtxAZ,q12,q22,q32,q1q2,
     -           debug,1)
            call az_vtxSUSY(a1vtxAZ2,a2vtxAZ2,a3vtxAZ2,q22,q12,q32,q1q2,
     -           debug,2)    

*            call zh_prpSUSY(a1zh,q12,q22,q32,q1q2,debug,HiggsType)
*            call gamh_prpSUSY(a1gamh,q12,q22,q32,q1q2,debug,HiggsType)

            if (HiggsType .ne. 3) then
               call ww_selfSUSY(a1slfW,q12,q22,q32,q1q2,debug)       
               call zz_selfSUSY(a1slfZ,q12,q22,q32,q1q2,debug)  
               call az_selfSUSY(a1slfAZ,q12,q22,q32,q1q2,debug,1) 
               call az_selfSUSY(a1slfAZ2,q22,q12,q32,q1q2,debug,2)
            end if
         end if

      end if  


** Resetting to corrected Higgs masses
      if ((model .eq. 2) .and. (MHloops .eq. 0)) then
         do t = 1, 4
            MHiggs(t) = MHiggs_loop(t)
         end do
         
         Mh0 = MHiggs_loop(1)
         MHH = MHiggs_loop(2)
         MA0 = MHiggs_loop(3)
         MHp = MHiggs_loop(4)
         
         Mh02 = Mh0**2
         MHH2 = MHH**2
         MA02 = MA0**2
         MHp2 = MHp**2
      end if
** End of resetting Higgs masses


** Re-scaling the self energy corrections if the Higgs being considered is the
** heavy CP-even Higgs or if propagator corrections are applied at loop level
      if (model .eq. 2) then
         if (propLEVEL .ne. 2) then
            if (HiggsType .eq. 2) then
               a1slfW = a1slfW*CBA/SBA
               a1slfZ = a1slfZ*CBA/SBA
               a1slfAZ = a1slfAZ*CBA/SBA
               a1slfAZ2 = a1slfAZ2*CBA/SBA
            end if
         else
            a1slfW = a1slfW*(Hgs(HiggsType,1) +
     -           Hgs(HiggsType,2)*CBA/SBA)
            a1slfZ = a1slfZ*(Hgs(HiggsType,1) +
     -           Hgs(HiggsType,2)*CBA/SBA)
            a1slfAZ = a1slfAZ*(Hgs(HiggsType,1) + 
     -           Hgs(HiggsType,2)*CBA/SBA)
            a1slfAZ2 = a1slfAZ2*(Hgs(HiggsType,1) +
     -           Hgs(HiggsType,2)*CBA/SBA)
         end if
      end if
* for fermions:
      if (model .eq. 2) then
         if (propLEVEL .ne. 2) then
            if (HiggsType .eq. 2) then
               a1selfWfm = a1selfWfm*CBA/SBA
               a1selfZfm = a1selfZfm*CBA/SBA
               a1selfAZfm = a1selfAZfm*CBA/SBA
               a1selfAZ2fm = a1selfAZ2fm*CBA/SBA
            end if
         else
            a1selfWfm = a1selfWfm*(Hgs(HiggsType,1) +
     -           Hgs(HiggsType,2)*CBA/SBA)
            a1selfZfm = a1selfZfm*(Hgs(HiggsType,1) +
     -           Hgs(HiggsType,2)*CBA/SBA)
            a1selfAZfm = a1selfAZfm*(Hgs(HiggsType,1) + 
     -           Hgs(HiggsType,2)*CBA/SBA)
            a1selfAZ2fm = a1selfAZ2fm*(Hgs(HiggsType,1) +
     -           Hgs(HiggsType,2)*CBA/SBA)
         end if
      end if
 

** Adding up all contributions to the form factors 
** ZZ fusion
      ahvvL(1,2,2) = a1vtxZ + a1slfZ 
      ahvvL(2,2,2) = a2vtxZ 
      ahvvL(3,2,2) = a3vtxZ 

** Z-photon fusion
      ahvvL(1,2,1) = a1vtxAZ2 + a1slfAZ2
*      ahvvL(1,2,1) = a1slfAZ2
      ahvvL(2,2,1) = a2vtxAZ2
      ahvvL(3,2,1) = a3vtxAZ2

** photon-Z fusion
      ahvvL(1,1,2) = a1vtxAZ + a1slfAZ
*      ahvvL(1,1,2) = a1slfAZ
      ahvvL(2,1,2) = a2vtxAZ
      ahvvL(3,1,2) = a3vtxAZ

** photon-photon fusion
      ahvvL(1,1,1) = a1vtxA
      ahvvL(2,1,1) = a2vtxA
      ahvvL(3,1,1) = a3vtxA

** WW fusion      
      ahvvL(1,3,4) = a1vtxW + a1slfW + a1zh + a1gamh
      ahvvL(2,3,4) = a2vtxW
      ahvvL(3,3,4) = a3vtxW

** If Higgs propagator corrections are included as an additional loop 
** correction
      if (propLevel .eq. 3) then
         ahvvL(1,2,2) = ahvvL(1,2,2) + a1propZ
         ahvvL(1,3,4) = ahvvL(1,3,4) + a1propW
      end if

      if (ewNLO .eq. 0) then
         ahvvL(1,2,2) = treeZ(HiggsType)
         ahvvL(1,3,4) = treeW(HiggsType)
      end if
         

      ahvvL(1,4,3) = ahvvL(1,3,4)
      ahvvL(2,4,3) = ahvvL(2,3,4)
      ahvvL(3,4,3) = ahvvL(3,3,4)


* for FERMIONS:** Adding up all contributions to the form factors 
** ZZ fusion
      ahvvFM(1,2,2) = a1vtxZfm + a1selfZfm
      ahvvFM(2,2,2) = a2vtxZfm 
      ahvvFM(3,2,2) = a3vtxZfm 
** Z-photon fusion
      ahvvFM(1,2,1) = a1vtxAZH2fm + a1selfAZ2fm
      ahvvFM(2,2,1) = a2vtxAZH2fm
      ahvvFM(3,2,1) = a3vtxAZH2fm
** photon-Z fusion
      ahvvFM(1,1,2) = a1vtxAZHfm + a1selfAZfm
      ahvvFM(2,1,2) = a2vtxAZHfm
      ahvvFM(3,1,2) = a3vtxAZHfm
** photon-photon fusion
      ahvvFM(1,1,1) = a1vtxAHfm
      ahvvFM(2,1,1) = a2vtxAHfm
      ahvvFM(3,1,1) = a3vtxAHfm
** WW fusion      
      ahvvFM(1,3,4) = a1vtxWfm + a1selfWfm 
      ahvvFM(2,3,4) = a2vtxWfm
      ahvvFM(3,3,4) = a3vtxWfm
** 'opposite' WW      
      ahvvFM(1,4,3) = ahvvFM(1,3,4)
      ahvvFM(2,4,3) = ahvvFM(2,3,4)
      ahvvFM(3,4,3) = ahvvFM(3,3,4)
* end of FOR FERMIONS


** Determining whether to include |M_{loop}|^2
      if (loopsqr_switchIN .and. (model .eq. 2)) then
         if ((treeZ(HiggsType) + a1propZ) .eq. 0D0) then
            loopSqr_switch = .true.
         else if ((treeW(HiggsType) + a1propW) .eq. 0D0) then
            loopSqr_switch = .true.
         else if (abs(ahvvFM(1,2,2)/(treeZ(HiggsType) +
     -           a1propZ)) .gt. 0.1D0) then
            loopSqr_switch = .true.
         else if (abs(ahvvFM(1,3,4)/(treeW(HiggsType) +
     -           a1propW)) .gt. 0.1D0) then
            loopSqr_switch = .true.
         else
            loopSqr_switch = .false.
         end if
      else
         loopSqr_switch = .false.
      end if


** Debugging output and kill command
      if ((debug .ne. 0) .and. (ewNLO .ne. 0)) then
*      if (ewNLO .ne. 0) then
         write(*,*)'  '

         write(*,*)'q12, q22, q32, q1q2 =', q12, q22, q32, q1q2 

         if (HiggsType .ne. 3) then
            write(*,*)'bornW, Z =', treeW(HiggsType) + a1propW, 
     -           treeZ(HiggsType) + a1propZ
         else
            write(*,*)'bornW, Z =', a1propW, a1propZ
         end if

         write(*,*)'a1propW, Z =', a1propW, a1propZ

         write(*,*)'a1selfW, Z, AZ =', a1slfW, a1slfZ, a1slfAZ

         write(*,*)'a1vtxW, Z =', a1vtxW, a1vtxZ
         write(*,*)'a1vtxA, AZ =', a1vtxA, a1vtxAZ

         write(*,*)'a1zh, a1gamh =', a1zh, a1gamh
         
         write(*,*)'a1factor - ZZ current =', ahvvL(1,2,2)
         write(*,*)'a2factor - ZZ current =', ahvvL(2,2,2)
         write(*,*)'a3factor - ZZ current =', ahvvL(3,2,2)

         write(*,*)'a1factor - ZA current =', ahvvL(1,2,1)
         write(*,*)'a2factor - ZA current =', ahvvL(2,2,1)
         write(*,*)'a3factor - ZA current =', ahvvL(3,2,1)

         write(*,*)'a1factor - AZ current =', ahvvL(1,1,2)
         write(*,*)'a2factor - AZ current =', ahvvL(2,1,2)
         write(*,*)'a3factor - AZ current =', ahvvL(3,1,2)

         write(*,*)'a1factor - AA current =', ahvvL(1,1,1)
         write(*,*)'a2factor - AA current =', ahvvL(2,1,1)
         write(*,*)'a3factor - AA current =', ahvvL(3,1,1)

         write(*,*)'a1factor - charged current =', ahvvL(1,3,4)
         write(*,*)'a2factor - charged current =', ahvvL(2,3,4)
         write(*,*)'a3factor - charged current =', ahvvL(3,3,4)


         write(*,*)'AND FOR FERMIONS:'
         write(*,*)'a1factor - ZZ current =', ahvvFM(1,2,2)
         write(*,*)'a2factor - ZZ current =', ahvvFM(2,2,2)
         write(*,*)'a3factor - ZZ current =', ahvvFM(3,2,2)

         write(*,*)'a1factor - ZA current =', ahvvFM(1,2,1)
         write(*,*)'a2factor - ZA current =', ahvvFM(2,2,1)
         write(*,*)'a3factor - ZA current =', ahvvFM(3,2,1)

         write(*,*)'a1factor - AZ current =', ahvvFM(1,1,2)
         write(*,*)'a2factor - AZ current =', ahvvFM(2,1,2)
         write(*,*)'a3factor - AZ current =', ahvvFM(3,1,2)

         write(*,*)'a1factor - AA current =', ahvvFM(1,1,1)
         write(*,*)'a2factor - AA current =', ahvvFM(2,1,1)
         write(*,*)'a3factor - AA current =', ahvvFM(3,1,1)

         write(*,*)'a1factor - charged current =', ahvvFM(1,3,4)
         write(*,*)'a2factor - charged current =', ahvvFM(2,3,4)
         write(*,*)'a3factor - charged current =', ahvvFM(3,3,4)

         write(*,*)'loopSqr_switch =', loopSqr_switch
         
*         if (sector .ne. 3) stop            
         stop
      end if
      

      end


********************************************************************************
********************************************************************************

** The following subroutines are all concerned with the calculation of the 
** matrix elements squared: i.e. the box and pentagon diagrams

********************************************************************************
********************************************************************************

*** This subroutine calculates the squared matrix elements of the boxes and 
*** pentagons.  It calls both the charged and neutral currents.

      subroutine boxes(pbar,fsign,uucc,uuss,ddcc,ddss,udsc,ducs,tree,
     &     fac,debug)

      implicit none


      double precision tree(6)


** i is the particle label according to the conventions
**                 i = 1         upper incoming quark
**                   = 2         upper outgoing quark
**                   = 3         lower incoming quark
**                   = 4         lower outgoing quark
**                   = 5         higgs
** mu is the spacetime index 0 -> 3.
      integer i, mu, j

** cur is the quark current label, where
**               cur = 1         uucc
**                   = 2         uuss
**                   = 3         ddcc
**                   = 4         ddss
**                   = 5         udsc
**                   = 6         ducs
** NC and CC are the neutral current (cur = 1->4) and charged current (cur = 
** 5->6) labels respectively.  See notebook 7, page 2.
      integer cur, NC, CC


** pbar and fsign are the inputs from the main vbfnlo program - the momentum
** and the particle/antiparticle label.  The final entry in each array is
** particle number i.  p is the physical momentum.  For pbar, i does not
** necessarily correspond to the particle label above - see notebook 7, page 12
      double precision pbar(0:3,6), p(0:3,5), intP(0:3,4)
      double precision intPmat(4,4), det
      integer fsign(6)


** massV contains the particle masses of the particular current that is being
** considered, and is labelled by particle number.  See notebook 7, page 1 for
** notation.
      double precision massV(5)


** qk is the quark/antiquark pattern label. See notebook 7, page 1.  It is 
** dependent on fsign: 
**               qk = 1          quark-quark
**                  = 2          antiquark-quark
**                  = 3          quark-antiquark
**                  = 4          antiquark-antiquark
      integer qk


** The momentum invariants
#include "global.inc"
      double precision SInv(n_p+n_v,n_p+n_v), TInv(n_p+n_v,n_p+n_v)


** The squared matrix elements.  The first entry in the array is 'cur' and the
** second is tree level (1) and loop level (2).  uucc etc. are the tree + loop
** values for each current's squared matrix element.
      double precision SqMatEl(6,2)
      double precision uucc, uuss, ddcc, ddss, udsc, ducs


** This is the multiplicative factor for the matrix elements.  It is calculated
** in qqhqq.F.
      double precision fac, facUSED
      double complex propFAC
      logical matFAC

*      double precision mu154, mu254, mu354, mu153, mu253, mu353
*      double precision mu143, mu243, mu343
*      integer meMH, meMH2
      integer incon


** External quark mass parameter - switch to 1 if you want the external quarks
** to have masses, thus allowing comparison (and debugging) with FormCalc, or 0
** if you want no external quark masses
      integer extQK
      parameter (extQK = 0)


** Do loop parameter
      integer type

** Debugging flag
      integer debug


** Declarations of the parameters
#include "mssm.inc"


** Uncomment for debugging output
*      write(*,*)'   '
*      write(*,*)'pre boxes...'
*      write(*,*)'uucc =', uucc
*      write(*,*)'uuss =', uuss
*      write(*,*)'ddcc =', ddcc
*      write(*,*)'ddss =', ddss
*      write(*,*)'udsc =', udsc
*      write(*,*)'ducs =', ducs


      incon = 0


** Setting to tree level Higgs masses and switching for sin-theta-W conventions
      if (model .eq. 2) then
         do type = 1, 4
            MHiggs(type) = MHiggs_tree(type)
         end do
         
         Mh0 = MHiggs(1)
         MHH = MHiggs(2)
         MA0 = MHiggs(3)
         MHp = MHiggs(4)

         Mh02 = Mh0**2
         MHH2 = MHH**2
         MA02 = MA0**2
         MHp2 = MHp**2

         SW = -SW
      end if
** End of setting tree level Higgs masses


** Inputs - the masses.  See notebook 7, page 1 for notation.
      do i = 1, 4
         massV(i) = 0D0
         massV(5) = MH
      end do
      
      
** Inputs - finding the quark pattern from fsign.  See notebook 7, page 1.
      if (fsign(1) .gt. 0) then
         if (fsign(3) .gt. 0) then
            qk = 1
         else
            qk = 3
         end if 
      else
         if (fsign(3) .gt. 0) then
            qk = 2
         else
            qk = 4
         end if
      end if


** Setting the physical momentum.  I realise the momentum notation is
** confusing - see Notebook 7, page 12.
      do i = 1, 4
         do mu = 0, 3
            p(mu,i) = pbar(mu,i)
         end do
      end do

      if (fsign(1) .eq. -1) then
         do mu = 0,3
            p(mu,1) = pbar(mu,2)
            p(mu,2) = pbar(mu,1)
         end do
      end if      
      if (fsign(3) .eq. -1) then
         do mu = 0,3
            p(mu,3) = pbar(mu,4)
            p(mu,4) = pbar(mu,3)
         end do
      end if

      do mu = 0, 3
         p(mu,5) = p(mu,1) + p(mu,3) - p(mu,2) - p(mu,4)
      end do


** Debugging output
      if (debug .ne. 0) then
         write(*,*)'  '
         do i = 1, 5
            do mu = 0, 3
               write(*,*)'i, mu, p(mu,i),pbar(mu,i) =', i, mu, p(mu,i), 
     -              pbar(mu,i)
            end do
         end do
      end if

c$$$      do mu = 0, 3
c$$$         intP(mu,1) = p(mu,3)
c$$$         intP(mu,2) = p(mu,3) - p(mu,4)
c$$$         intP(mu,3) = p(mu,3) - p(mu,4) - p(mu,5)
c$$$         intP(mu,4) = -p(mu,1)
c$$$      end do
c$$$
c$$$      do i = 1, 4
c$$$         do j = 1, 4
c$$$            intPmat(i,j) = 0D0
c$$$         end do
c$$$      end do
c$$$
c$$$      do i = 1, 4
c$$$         do j = 1, 4
c$$$            intPmat(i,j) = 2D0*intP(0,i)*intP(0,j)
c$$$            do mu = 1, 3
c$$$               intPmat(i,j) = intPmat(i,j) - (2D0*intP(mu,i)*intP(mu,j))
c$$$            end do
c$$$         end do
c$$$      end do
c$$$
c$$$
c$$$      det =  intPmat(1,1)*(intPmat(2,2)*(intPmat(3,3)*intPmat(4,4)-
c$$$     &     intPmat(3,4)*intPmat(4,3))+intPmat(2,3)*(intPmat(3,4)*
c$$$     &     intPmat(4,2)-intPmat(3,2)*intPmat(4,4))+intPmat(2,4)*
c$$$     &     (intPmat(3,2)*intPmat(4,3)-intPmat(3,3)*intPmat(4,2)))-
c$$$     &     intPmat(1,2)*(intPmat(2,1)*(intPmat(3,3)*intPmat(4,4)-
c$$$     &     intPmat(3,4)*intPmat(4,3))+intPmat(2,3)*(intPmat(3,4)*
c$$$     &     intPmat(4,1)-intPmat(3,1)*intPmat(4,4))+intPmat(2,4)*
c$$$     &     (intPmat(3,1)*intPmat(4,3)-intPmat(3,3)*intPmat(4,1)))+
c$$$     &     intPmat(1,3)*(intPmat(2,1)*(intPmat(3,2)*intPmat(4,4)-i
c$$$     &     ntPmat(3,4)*intPmat(4,2))+intPmat(2,2)*(intPmat(3,4)*
c$$$     &     intPmat(4,1)-intPmat(3,1)*intPmat(4,4))+intPmat(2,4)*
c$$$     &     (intPmat(3,1)*intPmat(4,2)-intPmat(3,2)*intPmat(4,1)))-
c$$$     &     intPmat(1,4)*(intPmat(2,1)*(intPmat(3,2)*intPmat(4,3)-
c$$$     &     intPmat(3,3)*intPmat(4,2))+intPmat(2,2)*(intPmat(3,3)*
c$$$     &     intPmat(4,1)-intPmat(3,1)*intPmat(4,3))+intPmat(2,3)*
c$$$     &     (intPmat(3,1)*intPmat(4,2)-intPmat(3,2)*intPmat(4,1)))
c$$$
c$$$      if (det**2 .lt. 10000D0) then
c$$$         write(*,*)'  '
c$$$         write(*,*)'det = ', det
c$$$         stop
c$$$      end if


** Calculating the momentum invariants and spinors.  This can be run OUTSIDE  
** the current do loop IFF the external quark masses are set to zero.
      if (extQK .eq. 0) then
         call momentum(cur, massV, p, SInv, TInv, debug)      
      end if

      if (debug .gt. 0) then
         write(*,*)'qk =', qk      
      end if


** Multiplying by SUSY factor when working in the MSSM: this consists of
** the (e.g.) sin(beta-alpha) difference in the tree level couplings and the
** propagator corrections if they are included at LO
      propFAC = (1D0,0D0)
      if (model .eq. 1) then
         facUSED = fac
         matFAC = .false.
      else if (model .eq. 2) then
         if (propLEVEL .eq. 1) then
            propFAC = (1D0,0D0)
            matFAC = .true.
         else if (propLEVEL .eq. 2) then
            propFAC = SBA*Hgs(HiggsType,1) + CBA*Hgs(HiggsType,2)
            matFAC = .false.
         else
            matFAC = .false.
            if (HiggsType .eq. 1) then
               propFAC = SBA
            else if (HiggsType .eq. 2) then
               propFAC = CBA
            end if
         end if
         facUSED = fac * propFAC * DCONJG(propFAC)
      end if

         
** Looping over the external quark currents
      do cur = 1, 6

         if (incon .eq. 0) then

** Setting the masses for the process 'cur'
            do i = 1, 5
               massV(i) = massCur(i,cur)
            end do

** Sorting out the charged current mass changes for different external quark 
** patterns.  See notebook 7, page 1 for details.
            if (cur .gt. 4) then
               if (qk .eq. 2 .or. qk .eq. 4) then
                  massV(1) = massCur(2,cur)
                  massV(2) = massCur(1,cur)
               end if
               if (qk .gt. 2) then
                  massV(3) = massCur(4,cur)
                  massV(4) = massCur(3,cur)
               end if
            end if


** Debugging output 
            if (debug .gt. 0) then
               do i = 1, 5
                  write(*,*)'i =', i
                  write(*,*)'massV =', massV(i)
               end do
            end if


** Calculating the momentum invariants and spinors.  This must be run inside  
** the current loop run IF the external quark masses are not set to zero.
            if (extQK .ne. 0) then
               call momentum(cur, massV, p, SInv, TInv, debug)
            end if  
           

** Clearing the squared matrix element values
            SqMatEl(cur,1) = 0D0
            SqMatEl(cur,2) = 0D0


** Neutral currents ...
            if (cur .le. 4) then
               NC = cur
** Calculating the squared matrix elements  
               call squaringNC(p,massV,qk,NC,SInv,TInv,SqMatEl,
     &              matFAC,debug)         
** Charged currents ...
            else if (cur .gt. 4) then
               CC = cur - 4
** Calculating the squared matrix elements
               call squaringCC(p,massV,qk,CC,SInv,TInv,SqMatEl,
     &              matFAC,debug)
** Charged currents ...
            end if

*            if (abs(((SqMatEl(cur,1)*facUSED)-tree(cur))*100/tree(cur)) 
*     &           .gt. 1D0) then
*               incon = 1
*               write(*,*)'tree(',cur,') =', tree(cur)
*               write(*,*)'SqMatEl(',cur,',1)*fac =', SqMatEl(cur,1)*fac
*            end if
*            if (abs(SqMatEl(cur,2)) .gt. 1000D0) then
*               incon = 1
*               write(*,*)'large loop'
*            end if

         end if                 ! end of incon if loop

** End of current loop
      end do


** Resetting to corrected Higgs masses and switching sin-theta-w conventions
      if (model .eq. 2) then
         do type = 1, 4
            MHiggs(type) = MHiggs_loop(type)
         end do
         
         Mh0 = MHiggs_loop(1)
         MHH = MHiggs_loop(2)
         MA0 = MHiggs_loop(3)
         MHp = MHiggs_loop(4)
         
         Mh02 = Mh0**2
         MHH2 = MHH**2
         MA02 = MA0**2
         MHp2 = MHp**2

         SW = -SW
      end if
** End of resetting Higgs masses


** serious debugging: outputs the previously calculated born level matrix
** element, the factor it's multiplied by, and the newly calculated born and
** loop contributions
*     write(*,*)'  '
*     write(*,*)'fac =', fac
*     write(*,*)'uucc =', tree(1)
*     write(*,*)'uuss =', tree(2)
*     write(*,*)'ddcc =', tree(3)
*     write(*,*)'ddss =', tree(4)
*     write(*,*)'udsc =', tree(5)
*     write(*,*)'ducs =', tree(6)

*     write(*,*)'SqMatEl(1,1) =', SqMatEl(1,1)*fac
*     write(*,*)'SqMatEl(2,1) =', SqMatEl(2,1)*fac
*     write(*,*)'SqMatEl(3,1) =', SqMatEl(3,1)*fac
*     write(*,*)'SqMatEl(4,1) =', SqMatEl(4,1)*fac
*     write(*,*)'SqMatEl(5,1) =', SqMatEl(5,1)*fac
*     write(*,*)'SqMatEl(6,1) =', SqMatEl(6,1)*fac

*      write(*,*)'SqMatEl(1,2) =', SqMatEl(1,2)*fac
*      write(*,*)'SqMatEl(2,2) =', SqMatEl(2,2)*fac
*      write(*,*)'SqMatEl(3,2) =', SqMatEl(3,2)*fac
*      write(*,*)'SqMatEl(4,2) =', SqMatEl(4,2)*fac
*      write(*,*)'SqMatEl(5,2) =', SqMatEl(5,2)*fac
*      write(*,*)'SqMatEl(6,2) =', SqMatEl(6,2)*fac
      

      if (incon .ne. 1) then
         uucc = uucc + (SqMatEl(1,2) * facUSED)
         uuss = uuss + (SqMatEl(2,2) * facUSED)
         ddcc = ddcc + (SqMatEl(3,2) * facUSED)
         ddss = ddss + (SqMatEl(4,2) * facUSED)
         udsc = udsc + (SqMatEl(5,2) * facUSED)
         ducs = ducs + (SqMatEl(6,2) * facUSED)
      end if


** Output for debugging: this will output the matrix squared element for each
** current and stop the code after one phase space point
      if (debug .gt. 3) then
         write(*,*)'qk =', qk
         write(*,*)'post boxes'
         write(*,*)'uucc =', uucc
         write(*,*)'uuss =', uuss
         write(*,*)'ddcc =', ddcc
         write(*,*)'ddss =', ddss
         write(*,*)'udsc =', udsc
         write(*,*)'ducs =', ducs
         write(*,*)'   '

         if (qk .eq. 4) then
            stop
         end if
      end if


      end


*******************************************************************************
*******************************************************************************
      
*** This subroutine calculates the invariants needed for the matrix elements
*** and calculates the spinor values so we can get the needed spinor products.
*** Most of this is stolen from FormCalc-5.2 (Thomas Hahn).

      subroutine momentum(cur, massV, p, SInv, TInv, debug)

      implicit none


** cur is the quark current label, where
**               cur = 1         uucc
**                   = 2         uuss
**                   = 3         ddcc
**                   = 4         ddss
**                   = 5         udsc
**                   = 6         ducs
** See notebook 7, page 2.
      integer cur


** Momentum and mass of the external particles input from the external 
** program.  Here:
**                    i = 1       initial upper quark
**                      = 2       final upper quark
**                      = 3       initial lower quark
**                      = 4       final lower quark
**                      = 5       higgs
** See notebook 7, page 1.
      double precision p(0:3,5), massV(5)


** Momentum and mass of the external particles.  q and massF are purely  
** internal definitions stolen from FormCalc).  q(t,i) and massF(i) where      
**                    i = 1       initial upper quark
**                      = 2       initial lower quark
**                      = 3       final upper quark
**                      = 4       higgs
**                      = 5       final lower quark
** and t = space-time index.  See notebook 7, page 1.
      double precision q(0:3,5), massF(5)


** Modulus of momentum of particle i and the unit vector describing its 
** direction
      double precision modq(5), vecq(3,5)


** The invariants SInv and TInv
#include "global.inc"
      double precision SInv(n_p+n_v,n_p+n_v), TInv(n_p+n_v,n_p+n_v)


** Do loop variables. i,j = particle number, mu = spacetime
      integer i, mu


** Debugging parameter
      integer debug


** Setting the internal definitions of momentum (i.e. switching to FormCalc
** notation from vbfnlo notation).
      do mu = 0, 3
         q(mu,1) = p(mu,1)
         q(mu,2) = p(mu,3)
         q(mu,3) = p(mu,2)
         q(mu,5) = p(mu,4)
         q(mu,4) = p(mu,5)         
      end do


** Setting the mass definitions (i.e. switching to FormCalc notation).  These
** are the masses of the EXTERNAL PARTICLES.
      massF(1) = massV(1)
      massF(2) = massV(3)
      massF(4) = massV(5)

      if (cur .le. 4) then
** i.e. neutral currents, particle mass in = particle mass out
         massF(3) = massF(1)
         massF(5) = massF(2)
      else
** i.e. charged currents, particle mass in .ne. particle mass out
         massF(3) = massV(2)
         massF(5) = massV(4)
      end if


** Setting masses of external quarks to zero
      do i = 1, 5
         massF(i) = 0D0
      end do
      massF(4) = massV(5)


** Debugging output
      if (debug .eq. 1) then
         do i = 1, 5
            write(*,*)'i =', i
            write(*,*)'ext.prog. massV(i) =', massV(i)
            write(*,*)'int. mom. massF(i) =', massF(i)
            do mu = 0, 3
               write(*,*)'i, mu =', i, mu
               write(*,*)'ext. prog. p(mu,i) =', p(mu,i)
               write(*,*)'int. mom. q(mu,i) =', q(mu,i)
            end do
         end do
      end if


** Calculating the modulus and direction of the momentum
      do i = 1, 5
         modq(i) = q(1,i)**2 + q(2,i)**2 + q(3,i)**2
         modq(i) = sqrt(modq(i))

         do mu = 1, 3
            vecq(mu,i) = q(mu,i) / modq(i)
         end do
      end do


** Now setting KK(i,x,y) and vec(x,y,k(i))  - see VecSet.F  
      do i = 1, 5
         call FC_kin(i,massF(i),modq(i),vecq(1,i),vecq(2,i),vecq(3,i))
      end do

      
** The invariants
      call FCinvariants(SInv,TInv)


** Debugging messages
      if (debug .eq. 2) then
         do i = 1, 5
            write(*,*)'     '
            write(*,*)'i = ', i

            do mu = 1, 3
               write(*,*)'mu, vecq(mu,i) =', mu, vecq(mu,i)
            end do

            write(*,*)'modq(i) =',  modq(i)
         end do
         stop
      endif

      
      end


*******************************************************************************
*******************************************************************************
*******************************************************************************
        
*** This subroutine calculates the subtraction term to deal with the real and
*** collinear divergences.  See notebook 7, page 14.  This uses the dipole
*** formalism, as described in hep-ph/9904440.

      subroutine realCor(mass,qk,cur,p,q,Dsub,debug)

      implicit none

    
** The subtraction term
      double precision Dsub


** Current label, where:
**          cur = 1            NC = 1 :     u + c -> u + c
**              = 2            NC = 2 :     u + s -> u + s
**              = 3            NC = 3 :     d + c -> d + c
**              = 4            NC = 4 :     d + s -> d + s
**              = 5            CC = 1 :     u + s -> d + c
**              = 6            CC = 2 :     d + c -> u + s
      integer cur

** qk is the quark/antiquark pattern label. It is dependent on fsign: 
**               qk = 1          quark-quark
**                  = 2          antiquark-quark
**                  = 3          quark-antiquark
**                  = 4          antiquark-antiquark
      integer qk 

** Momenta p(mu,i) and charge q(i) of the external quarks.  Here,
**          i = 1             incoming upper quark
**            = 2             outgoing upper quark
**            = 3             incoming lower quark
**            = 4             outgoing lower quark
      double precision p(0:3,5), q(4), chargeQ(4), sig(4)

** The momentum dot product p(i).p(j)
      double precision p2(4,4)      

** NOT USED in this version of the corrections
      double precision mass(4)


** Do loop variables for particle number
      integer i, j, k


** Logarithmic function used in the subtraction term (see notebook 7, page 14)
      double precision rlCfn
      external rlCfn


** Checking charge conservation
      double precision qCon

      integer debug


** Declaration of parameters
#include "mssm.inc"  
#include "global.inc"  


** Clearing the value of the subtraction term
      Dsub = 0D0


** Defining the chargeQ
      do i = 1,4
         chargeQ(i) = q(i)
      end do

** For neutral currents: upper chargeQ in = upper chargeQ out
      if (cur .le. 4) then
         chargeQ(2) = chargeQ(1)
         chargeQ(4) = chargeQ(3)
      end if


** Defining the quark / antiquarkness:  Sig is +1 for incoming fermion / 
** outgoing antifermion and sig is -1 for outgoing fermion / incoming 
** antifermion

** Qk = 1; upper incoming fermion
      sig(1) = 1D0
** Qk = 1; upper outgoing fermion
      sig(2) = -1D0
** Qk = 1; lower incoming fermion
      sig(3) = 1D0
** Qk = 1; upper outgoing fermion
      sig(4) = -1D0

      if (qk .eq. 2 .or. qk .eq. 4) then
         sig(1) = -sig(1)
         sig(2) = -sig(2)
      end if
      if (qk .eq. 3 .or. qk .eq. 4) then
         sig(3) = -sig(3)
         sig(4) = -sig(4)
      end if
      

** Checking charge conservation
      qCon = 0D0
      do i = 1, 4
         qCon = qCon + sig(i)*chargeQ(i)
      end do
      if (abs(qCon) .gt. 1D-5) then
         write(*,*)'Charge conservation is NOT WORKING'
         write(*,*)'sum of charges, subtraction term =', qCon
         stop
      end if


** Debugging output           
      if (debug .ne. 0) then
         write(*,*)'     '
         write(*,*)'cur, qk =', cur, qk
         do i = 1, 4
            write(*,*)'i =', i
            write(*,*)'chargeQ(i),q(i) =', chargeQ(i)
            write(*,*)'sig(i) =', sig(i)
            do j = 0, 3
               write(*,*)'j, p(j,i) =', j, p(j,i)
            end do
         end do
      end if
      

** Calculating the momentum dot product p2(i,j) = 2*p(i).p(j)
      do i = 1, 4
         do j = 1, 4
            p2(i,j) = (p(0,i)*p(0,j)) - (p(1,i)*p(1,j))  
            p2(i,j) = p2(i,j) - (p(2,i)*p(2,j))  
            p2(i,j) = p2(i,j) - (p(3,i)*p(3,j))
            p2(i,j) = 2D0*p2(i,j)
         end do
      end do 


** The subtraction term itself

      Dsub = chargeQ(2)*sig(2)*chargeQ(4)*sig(4)*
     -     ((2D0*rlCfn(p2(2,4)))-(2D0*(Pi**2)/3D0)+3D0)

      Dsub = Dsub + chargeQ(1)*sig(1)*chargeQ(3)*sig(3)*
     -     ((2D0*rlCfn(p2(1,3)))-(2D0*(Pi**2)/3D0)+4D0)

      Dsub = Dsub + chargeQ(1)*sig(1)*chargeQ(2)*sig(2)*
     -     ((2D0*rlCfn(p2(1,2)))-((Pi**2)/3D0)+0.5D0)

      Dsub = Dsub + chargeQ(3)*sig(3)*chargeQ(2)*sig(2)*
     -     ((2D0*rlCfn(p2(2,3)))-((Pi**2)/3D0)+0.5D0)

      Dsub = Dsub + chargeQ(1)*sig(1)*chargeQ(4)*sig(4)*
     -     ((2D0*rlCfn(p2(1,4)))-((Pi**2)/3D0)+0.5D0)
      
      Dsub = Dsub + chargeQ(3)*sig(3)*chargeQ(4)*sig(4)*
     -     ((2D0*rlCfn(p2(3,4)))-((Pi**2)/3D0)+0.5D0)

      Dsub = Dsub * (-AlfaQED/(2D0*Pi))


** Debugging output
      if (debug .ne. 0) then
         if (cur .eq. 1) then
            write(*,*)' '
         end if
         write(*,*)'qk, cur, Dsub =', qk, cur, Dsub
      end if


      end


*******************************************************************************
*******************************************************************************

*** This is the logarithmic function used when calculating the subtraction
*** term.  See notebook 7, page 14.

      double precision function rlCfn(p2)

      implicit none

** The momentum dot product
      double precision p2, r

** Declaration of parameters
#include "mssm.inc"


      r = p2

      if (r .lt. 0D0) then
         r = -r
      end if

      
      rlCfn = log((massQ**2)/r)*log((massG)/r) 
      
      rlCfn = rlCfn + log((massG)/r)    
      
      rlCfn = rlCfn - (0.5D0*((log((massQ**2)/r))**2))  
      
      rlCfn = rlCfn + (0.5D0*log((massQ**2)/r)) 
      

      end


*******************************************************************************


*******************************************************************************
*******************************************************************************
        
*** This subroutine calculates the subtraction term to deal with the IR 
*** divergences.  It uses a soft photon factor, as described in 
*** hep-ph/0709.1075, section 7.  It is stolen pretty much completely from
*** FormCalc.  
*** NOTE: This is not used by default, it's just here to allow an additional
*** check of the IR divergences.

      subroutine realCorZ(mass,qk,cur,mom,q,Dsub,debug)

      implicit none


** The subtraction term
      double precision Dsub


** Momenta p(mu,i) and chargeQ(i) of the external quarks.  To avoid mistakes,
** I use a completely internal definition for these: they are taken from the
** external mom and q. Here,
**          i = 1             incoming upper quark
**            = 2             outgoing upper quark
**            = 3             incoming lower quark
**            = 4             outgoing lower quark
      double precision mom(0:3,5), q(4), chargeQ(5), sig(4)

** The momentum dot product p(i).p(j)
      double precision p2(4,4), p(0:3,5)

** Mass of external particles:  mass is the external definition, m the internal
      double precision mass(5), m(5)


** Current label, where:
**          cur = 1            NC = 1 :     u + c -> u + c
**              = 2            NC = 2 :     u + s -> u + s
**              = 3            NC = 3 :     d + c -> d + c
**              = 4            NC = 4 :     d + s -> d + s
**              = 5            CC = 1 :     u + s -> d + c
**              = 6            CC = 2 :     d + c -> u + s
      integer cur

** qk is the quark/antiquark pattern label. It is dependent on fsign: 
**               qk = 1          quark-quark
**                  = 2          antiquark-quark
**                  = 3          quark-antiquark
**                  = 4          antiquark-antiquark
      integer qk 


      double precision PhaseSpaceIntegral
      external PhaseSpaceIntegral
      
** Do loop integers
      integer i, j
      
      double precision loglambda
      common /infrared/ loglambda
      
      double precision r

** Variable depending on centre of mass energy of collider, taken from ecm 
** which is set in input.dat)
      double precision ESOFTMAX

      integer debug

** Declaration of parameters
#include "mssm.inc"  
#include "global.inc"  


      ESOFTMAX = 0.1D0*ecm

** Setting the masses
      m(4) = mass(5)
      m(5) = mass(4)
      m(1) = mass(1)
      m(2) = mass(3)
      m(3) = mass(2)
        
** Defining the charge
      chargeQ(1) = q(1)
      chargeQ(2) = q(3)
      chargeQ(3) = q(2)
      chargeQ(4) = 0D0
      chargeQ(5) = q(4)

** Setting the momentum
      do i = 0, 3
         p(i,1) = mom(i,1)
         p(i,2) = mom(i,3)
         p(i,3) = mom(i,2)
         p(i,5) = mom(i,4)
         p(i,4) = mom(i,5)
      end do        


** Re-setting the charges and masses to the actual, physical values
      if (cur .le. 4) then
         chargeQ(3) = chargeQ(1)
         chargeQ(5) = chargeQ(2)
         m(3) = m(1)
         m(5) = m(2)
      end if  

      if (qk .eq. 2) then
         chargeQ(1) = -chargeQ(1)
         chargeQ(3) = -chargeQ(3)
      else if (qk .eq. 3) then
         chargeQ(2) = -chargeQ(2)
         chargeQ(5) = -chargeQ(5)
      else if (qk .eq. 4) then
         chargeQ(1) = -chargeQ(1)
         chargeQ(2) = -chargeQ(2)
         chargeQ(3) = -chargeQ(3)
         chargeQ(5) = -chargeQ(5)
      end if


** Debugging output
      if (debug .gt. 1) then
        write(*,*)'   '
        write(*,*)'real Cor ...'
        do i = 1, 5
           write(*,*)'i =', i
           write(*,*)'m(i) =', m(i)
           write(*,*)'chargeQ(i) =', chargeQ(i)
           write(*,*)'p(0,i), p(1,i) =', p(0,i), p(1,i)
           write(*,*)'p(2,i), p(3,i) =', p(2,i), p(3,i)
        end do        
      end if


      loglambda = log(4*(ESOFTMAX)**2/massG)


      Dsub = 0D0
      do i = 1, 5
         do j = i, 5
            r = chargeQ(i)*chargeQ(j)
            if( r .ne. 0 ) Dsub = Dsub +
     &           sign(1, 2 - i)*sign(1, 2 - j)*
     &           r*PhaseSpaceIntegral(m, p, i, j)
         enddo
      enddo
      Dsub = -4*pi*AlfaQED/(2*pi)**3*Dsub

** Debugging output
      if (debug .ne. 0) then
         if (cur .eq. 1) then
            write(*,*)'  '
         end if
         write(*,*)'Dsub =', Dsub
      end if

     
      end


*******************************************************************************


************************************************************************
** PhaseSpaceIntegral calculates the integral over the soft photon's
** phase space.  See A. Denner, Fortschr. d. Physik, 41 (1993) 4.
** NOTE: This is stolen point-blank from FormCalc

      double precision function PhaseSpaceIntegral(m, p, i, j)
      implicit none
      integer i, j

      double precision mi, ki0, ki, kidelta, kilog
      double precision mj, kj0, kj, kjdelta, kjlog
      double precision kikj, alpha, d, vi

      double precision loglambda
      common /infrared/ loglambda

      double precision m(5)
        
      double precision p(0:3,5)

      double complex SLi2
      external SLi2

#include "global.inc"


      mi = m(i)
      ki0 = p(0,i)
      ki = sqrt(p(1,i)**2 + p(2,i)**2 + p(3,i)**2)
      kidelta = (m(i)**2)/(ki0+ki)
      kilog = log(kidelta/(ki0 + ki))

      if( i .eq. j ) then
         if( ki .lt. 1D-14 ) then
	    PhaseSpaceIntegral = pi*(loglambda - 2)
         else
	    PhaseSpaceIntegral = pi*(loglambda + ki0/ki*kilog)
         endif
         return
      endif
      
      mj = m(j)
      kj0 = p(0,j)
      kj =  sqrt(p(1,j)**2 + p(2,j)**2 + p(3,j)**2)
      kjdelta = (m(j)**2)/(kj0+kj)
      kjlog = log(kjdelta/(kj0 + kj))

      kikj = p(0,i)*p(0,j) - p(1,i)*p(1,j) - p(2,i)*p(2,j) -
     &     p(3,i)*p(3,j)
      
      d = sqrt((kikj - mi*mj)*(kikj + mi*mj))
      alpha = (kikj + d)/mi**2
      if( alpha*ki0 .lt. kj0 ) alpha = (kikj - d)/mi**2
      
      d = 1/((alpha*mi - mj)*(alpha*mi + mj))
      vi = 2*(alpha*ki0 - kj0)*d
      
      PhaseSpaceIntegral = 4*pi*alpha*kikj*d * (
     &     log(mi/mj*abs(alpha))*loglambda +
     &     .25D0*(kilog - kjlog)*(kilog + kjlog) +
     &     SLi2(1 - (ki0 + ki)*alpha*vi) + SLi2(1 - kidelta*alpha*vi)-
     &     SLi2(1 - (kj0 + kj)*vi) - SLi2(1 - kjdelta*vi) )


	end


*******************************************************************************
*******************************************************************************

** This function is utterly stolen from FormCalc.

      double complex function SLi2(x)
      implicit none
      
      double precision x
        
      double precision pi12
      parameter (pi12 = .822467033424113218236207583323D0)

      double complex res, dummy
      integer ier, ipi12

      ier = 0
      call Sffzxdl(res, ipi12, dummy, x, -1, ier)

      SLi2 = res + ipi12*pi12

      end


*******************************************************************************
*******************************************************************************

** This subroutine is stolen completely from FormCalc.

*###[ ffzxdl:
      subroutine Sffzxdl(zxdilo,ipi12,zlog,x,ieps,ier)
***#[*comment:***************************************************
*	Computes the dilogarithm (Li2, Sp) for any (real) x	*
*	to precision precx. If an error message is given add	*
*	more bf's. For x > 1 the imaginary part is		*
*	 -/+i*pi*log(x), corresponding to x+ieps.		*
*	The number of factors pi^2/12 is passed separately in	*
*	ipi12 for accuracy.  We also calculate log(1-x)		*
*	which is likely to be needed.				*
*								*
*	Input:	x	(real)					*
*		ieps	(integer,+/-1)				*
*								*
*	Output: zxdilo	(complex) the dilog mod factors pi2/12	*
*		ipi12	(integer) these factors			*
*		zlog	(complex) log(1-x)			*
*								*
*	Calls:	log,Sdfflo1					*
*								*
***#]*comment:***************************************************
*  #[ declarations:

      implicit none

*	arguments
      integer ipi12,ieps,ier
      DOUBLE PRECISION x
      DOUBLE COMPLEX zxdilo,zlog

*	local variables
      integer jsgn
      DOUBLE PRECISION fact,u,u2,Sdfflo1,Sffbnd,a,xdilo,
     +     xprec,bdn02,bdn05,bdn10,bdn15
      DOUBLE COMPLEX cy,cfact
      external Sffbnd,Sdfflo1
      save xprec,bdn02,bdn05,bdn10,bdn15

#include "ff.h"


*  #] declarations:
*  #[ initialisations:
      data xprec /-1D0/
      if ( xprec .ne. precx ) then
         xprec = precx
         bdn02 = Sffbnd(1,2,bf)
         bdn05 = Sffbnd(1,5,bf)
         bdn10 = Sffbnd(1,10,bf)
         bdn15 = Sffbnd(1,15,bf)
      endif
*  #] initialisations:
*  #[ exceptional cases:
      if ( x .eq. 1) then
         zxdilo = 0
         zlog = -99999
         ipi12 = 2
         return
      elseif (x .eq. -1) then
         zxdilo = 0
         zlog = xlg2
         ipi12 = -1
         return
      elseif (x .eq. .5D0) then
         zxdilo = - xlg2**2/2
         zlog = -xlg2
         ipi12 = 1
         return
      elseif ( abs(x) .lt. precx ) then
         zxdilo = x
         zlog = -x
         ipi12 = 0
         return
      endif
*  #] exceptional cases:
*  #[ transform to (-1,.5):
      if (x .lt. -1) then
         fact = log(-x)
         cy = - fact**2/2
         ipi12 = -2
         if ( -x*xloss .gt. 1 ) then
            u = -Sdfflo1(1/x,ier)
         else
            u = -log(1-1/x)
         endif
         zlog = log(1-x)
         jsgn = -1
      elseif ( x .lt. .5D0) then
         cy = 0
         ipi12 = 0
         if ( abs(x) .lt. xloss ) then
            zlog = Sdfflo1(x,ier)
         else
            zlog = log(1-x)
         endif
         u = -DBLE(zlog)
         jsgn = 1
      elseif ( x .le. 2 ) then
         u = -log(x)
         if ( abs(1-x) .lt. xalogm ) then
            cy = 0
         elseif ( x .lt. 1 ) then
            zlog = log(1-x)
            cy = DBLE(u)*zlog
         elseif ( ieps .gt. 0 ) then
            zlog = DCMPLX(log(x-1),-pi)
            cy = DBLE(u)*zlog
         else
            zlog = DCMPLX(log(x-1),+pi)
            cy = DBLE(u)*zlog
         endif
         ipi12 = 2
         jsgn = -1
      else
         if ( ieps .gt. 0 ) then
            cfact = DCMPLX(log(x),-pi)
            zlog = DCMPLX(log(x-1),-pi)
         else
            cfact = DCMPLX(log(x),+pi)
            zlog = DCMPLX(log(x-1),+pi)
         endif
         cy = - cfact**2/2
         ipi12 = -2
         if ( x*xloss .gt. 1 ) then
            u = -Sdfflo1(1/x,ier)
         else
            u = -log(1-1/x)
         endif
         jsgn = -1
      endif
*  #] transform to (-1,.5):
*  #[ calculate dilog:
      if ( abs(u) .lt. xalog2 ) then
         xdilo = u
      else
         u2 = u**2
         a = abs(u2)
         if ( a .gt. bdn15 ) then
	    xdilo = u2*(bf(16) + u2*(bf(17) + u2*(bf(18) +
     +           u2*(bf(19) + u2*bf(20) ))))
         else
	    xdilo = 0
         endif
         if ( a .gt. bdn10 ) then
	    xdilo = u2*(bf(11) + u2*(bf(12) + u2*(bf(13) +
     +           u2*(bf(14) + u2*(bf(15) + xdilo)))))
         endif
         if ( a .gt. bdn05 ) then
	    xdilo = u2*(bf(6) + u2*(bf(7) + u2*(bf(8) +
     +           u2*(bf(9) + u2*(bf(10) + xdilo)))))
         endif
         if ( a .gt. bdn02 ) then
	    xdilo = u2*(bf(3) + u2*(bf(4) + u2*(bf(5) + xdilo)))
         endif
*     watch the powers of u.
         xdilo = u + u2*(bf(1) + u*(bf(2) + xdilo))
      endif
      if(jsgn.eq.1)then
         zxdilo =  DBLE(xdilo) + cy
      else
         zxdilo = -DBLE(xdilo) + cy
      endif
*  #] calculate dilog:
*###] ffzxdl:
      end


********************************************************************************


*###[ Sffbnd:
	DOUBLE PRECISION function Sffbnd(n1,n2,array)
*************************************************************************
*									*
*	calculate bound = (precx*|a(n1)/a(n1+n2)|^(1/n2) which is the	*
*	maximum value of x in a series expansion sum_(i=n1)^(n1+n2)	*
*	a(i)*x(i) to give a result of accuracy precx (actually of |next	*
*	term| < prec							*
*									*
*************************************************************************
	implicit none
	integer n1,n2
	DOUBLE PRECISION array(n1+n2)
	include 'ff.h'
	if ( array(n1+n2) .eq. 0 ) then
	   print *,'Sffbnd: fatal: array not initialized; did you call ',
     +		'ffini?'
	   stop
	endif
	Sffbnd = (precx*abs(array(n1)/array(n1+n2)))**(1/DBLE(n2))
*###] Sffbnd:
	end


********************************************************************************


*###[ Sdfflo1:
	DOUBLE PRECISION function Sdfflo1(x,ier)
***#[*comment:***************************************************
*	calculates log(1-x) for |x|<.14 in a faster way to ~15	*
*	significant figures.					*
***#]*comment:***************************************************
*  #[ declarations:
	implicit none
	integer ier
	DOUBLE PRECISION x,bdn01,bdn05,bdn10,bdn15,bdn19,xprec,
     +		xa,Sffbnd
	DOUBLE COMPLEX Szxfflg
	external Sffbnd,Szxfflg
	save xprec,bdn01,bdn05,bdn10,bdn15,bdn19
	include 'ff.h'
*  #] declarations:
*  #[ initialisation:
	data xprec /-1D0/
	if ( xprec .ne. precx ) then
	    xprec = precx
*	    determine the boundaries for 1,5,10,15 terms
	    bdn01 = Sffbnd(1,1,xninv)
	    bdn05 = Sffbnd(1,5,xninv)
	    bdn10 = Sffbnd(1,10,xninv)
	    bdn15 = Sffbnd(1,15,xninv)
	    bdn19 = Sffbnd(1,19,xninv)
	endif
*  #] initialisation:
*  #[ calculations:
	xa = abs(x)
	if ( xa .gt. bdn19 ) then
	    Sdfflo1 = DBLE(Szxfflg(1-x,0,0D0,ier))
	    return
	endif
	if ( xa .gt. bdn15 ) then
	    Sdfflo1 = x*( xninv(16) + x*( xninv(17) + x*( xninv(18) +
     +		x*( xninv(19) + x*xninv(20) ))))
	else
	    Sdfflo1 = 0
	endif
	if ( xa .gt. bdn10 ) then
	    Sdfflo1 = x*( xninv(11) + x*( xninv(12) + x*( xninv(13) +
     +		x*( xninv(14) + x*( xninv(15) + Sdfflo1 )))))
	endif
	if ( xa .gt. bdn05 ) then
	    Sdfflo1 = x*( xninv(6) + x*( xninv(7) + x*( xninv(8) +
     +		x*( xninv(9) + x*( xninv(10) + Sdfflo1 )))))
	endif
	if ( xa .gt. bdn01 ) then
	    Sdfflo1 = x*( xninv(2) + x*( xninv(3) + x*( xninv(4) +
     +		x*( xninv(5) + Sdfflo1 ))))
	endif
	Sdfflo1 = - x*( xninv(1) + Sdfflo1 )
*  #] calculations:
*###] Sdfflo1:
	end


********************************************************************************


*###[ Szxfflg:
	DOUBLE COMPLEX function Szxfflg(x,ieps,y,ier)
***#[*comment:***********************************************************
*									*
*	Calculate the complex logarithm of x.  The following cases	*
*	are treted separately:						*
*		|x| too small:		give warning and return 0	*
*					(for Absoft, Apollo DN300)	*
*		|x| < 0:		take sign according to ieps	*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
*
*	arguments
*
	implicit none
	integer ieps,ier
	DOUBLE PRECISION x,y
*
*	local variables
*
	DOUBLE PRECISION xlog
*
*	common blocks
*
	include 'ff.h'
*  #] declarations:
*  #[ calculations:
	if ( abs(x) .lt. xalogm ) then
	    Szxfflg = 0
	elseif ( x .gt. 0 ) then
	    Szxfflg = log(x)
	else
	    xlog = log(-x)
*	    checked imaginary parts 19-May-1988
	    if ( abs(ieps) .eq. 1 ) then
		if ( y*ieps .lt. 0 ) then
		    Szxfflg = DCMPLX(xlog,-pi)
		else
		    Szxfflg = DCMPLX(xlog,pi)
		endif
	    elseif ( ieps .eq. 2 ) then
		Szxfflg = DCMPLX(xlog,-pi)
	    elseif ( ieps .eq. -2 ) then
		Szxfflg = DCMPLX(xlog,+pi)
	    else
*		call fferr(52,ier)
		Szxfflg = DCMPLX(xlog,pi)
	    endif
	endif
*  #] calculations:
*###] Szxfflg:
	end


********************************************************************************
