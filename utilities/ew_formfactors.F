*******************************************************************************
*******************************************************************************

******************************************************************
***  ew_formfactors.F                                          ***
***  sophy@particle.uni-karlsruhe.de                           ***
***                                                            ***
***  This set of subroutines contains the contributions to     ***
***  the formfactors, for the STANDARD MODEL and SUSY.         ***
******************************************************************

*******************************************************************************
*******************************************************************************

*** This subroutine calculates the tree-level part of the a1 form factor for 
*** the T-mu-nu anomalous coupling for the VV-higgs vertex in the STANDARD 
*** MODEL

      subroutine treeSM(treeW,treeZ,debug,higgno)

      implicit none


** Value of tree level a1 coupling for W and Z bosons coupling.  The first 
** entry in each array is the light higgs, the second entry is the heavy higgs.
      double complex treeW(3), treeZ(3)


** higgno defines which of the neutral higgs particles we are considering.  
** In the SM the mass of the higgs can be set to match the mass of a higgs in 
** the defined susy scenario. 
      integer higgno

      integer debug


** Parameter declaration
#include "mssm.inc"

      treeW(higgno) = MW*EL/SW
      treeZ(higgno) = 1/CW2*(MW*EL/SW) 


      if (debug .ne. 0) then
         write(*,*)'WWH tree=', treeW(higgno)
         write(*,*)'ZZH tree=', treeZ(higgno)
      end if

      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the WWH vertex in 
*** the standard model 

      subroutine ww_vtxSM(a1vtxW,a2vtxW,a3vtxW,q12,q22,q32,q1q2,
     -     debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** The arguments of the integrals
      integer*8 iint1(3), iint2(3), iint3(3), iint4(3)
      integer*8 iintL1(3), iintL2(3)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6, iintB7 
      integer*8 iintB8, iintB9, iintB10, iintB11, iintB12, iintB13
      integer*8 iintB14


** These are the contributions from the counter-term (g-mu-nu) and the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1, loopa2, loopa3

** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for WW-higgs 
      double complex a1vtxW, a2vtxW, a3vtxW


** Particle do loop parameters 
      integer g4ini, Gen4

** em coupling constant
      double precision Alfa, Alfa2


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

      integer debug


** This is needed for looptools 
#include "looptools.h"

** Parameter declarations
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 quarks, 2 for all fermions and 3 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


** Now for the arguments of the integrals - first the quarks
      do Gen4 = g4ini, 3
         iint1(Gen4) = Bget(q12,Mf2(3,Gen4),Mf2(4,Gen4))
         iint2(Gen4) = Cget(q32,q12,q22,Mf2(3,Gen4),
     -        Mf2(3,Gen4),Mf2(4,Gen4))       
         iint3(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))         
         iint4(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),
     -        Mf2(4,Gen4),Mf2(4,Gen4))
      enddo

** The arguments of the leptonic integrals
      if (sector .ne. 1) then
         do Gen4 = 1,3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))            
            iintL2(Gen4) = Cget(q22,q32,q12,0.D0,Mf2(2,Gen4),Mf2(2,Gen4))
         enddo
      end if

** The arguments of the bosonic integrals
      if (sector .eq. 3) then
         iintB1 = Bget(q32,MH2,MH2)         
         iintB2 = Bget(q32,MW2,MW2)
         iintB3 = Bget(q32,MZ2,MZ2)
         iintB4 = Bget(q12,MW2,MZ2)
         iintB5 = Bget(q12,0.D0,MW2)
         iintB6 = Bget(q22,0.D0,MW2)
         iintB7 = Bget(q12,MH2,MW2)
         iintB8 = Bget(q22,MH2,MW2)
         iintB9 = Bget(q22,MW2,MZ2)
         iintB10 = Cget(q32,q12,q22,MH2,MH2,MW2)
         iintB11 = Cget(q32,q12,q22,MW2,MW2,MZ2)
         iintB12 = Cget(q22,q32,q12,0.D0,MW2,MW2)
         iintB13 = Cget(q22,q32,q12,MH2,MW2,MW2)
         iintB14 = Cget(q22,q32,q12,MW2,MZ2,MZ2)
      end if


** Reset each component to zero, just in case
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0)
      ct = (0D0, 0D0)


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
      ct = MW/(2.D0*MW2*SW2)*(EL*(SW2*dMWsq1 + 2*MW2*SW2*dZet1 + 
     -     MW2*SW2*dZH1 + 2*MW2*SW2*dZW1 - 2*MW2*(dSW1*SW)))/SW

** counterterm for the fermions
      ctFM = MW/(2.D0*MW2*SW2)*(EL*(SW2*dMWsq1fm + 2*MW2*SW2*dZet1fm + 
     -     MW2*SW2*dZH1fm + 2*MW2*SW2*dZW1fm - 2*MW2*(dSW1fm*SW)))/SW
      
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'WWH ct =', ct
         write(*,*)'WWH fermion ct =', ctFM
      end if


** Now the loop
      
** Quark contributions
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + 1/(16.D0*MW*Pi*SW2)*(Cval(cc0,iint2(Gen4))*
     -        (-(3*Alfa*((q12-q32-q22)*EL*Mf2(3,Gen4))) + 
     -        12*Alfa*(EL*Mf2(3,Gen4)**2)))/SW + 
     -        (3*Alfa)/(4.D0*MW*Pi*SW2)*
     -        (EL*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4))/SW - 
     -        (3*Alfa)/(2.D0*MW*Pi*SW2)*
     -        (EL*Cval(cc00,iint2(Gen4))*Mf2(3,Gen4))/SW - 
     -        (3*Alfa)/(8.D0*MW*Pi*SW2)*((q12-q22-2*q32)*EL*
     -        Cval(cc1,iint2(Gen4))*Mf2(3,Gen4))/SW - 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*((q12-q32-5*q22)*EL*
     -        Cval(cc2,iint2(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa)/(4.D0*MW*Pi*SW2)*
     -        (EL*Bval(bb0,iint3(Gen4))*Mf2(4,Gen4))/SW - 
     -        (3*Alfa)/(2.D0*MW*Pi*SW2)*
     -        (EL*Cval(cc00,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*((q12+3*q22-q32)*
     -        EL*Cval(cc1,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (9*Alfa)/(16.D0*MW*Pi*SW2)*((q12-1/3.D0*q32+1/3.D0*q22)*
     -        EL*Cval(cc2,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (3*Alfa)/(4.D0*MW*Pi*SW2)*
     -        (EL*Cval(cc0,iint4(Gen4))*Mf2(3,Gen4)*Mf2(4,Gen4))/SW

         loopa2 = loopa2 + 1/(16.D0*MW*Pi*SW2)*(Cval(cc0,iint2(Gen4))*
     -        (-(3*Alfa*(2*EL*Mf2(3,Gen4))) + 
     -        12*Alfa*(EL*0*Mf2(3,Gen4)**2)))/SW - 
     -        (3*Alfa)/(8.D0*MW*Pi*SW2)*
     -        (4*EL*Cval(cc1,iint2(Gen4))*Mf2(3,Gen4))/SW - 
     -        (3*Alfa)/(2.D0*MW*Pi*SW2)*
     -        (EL*(Cval(cc11,iint2(Gen4)) + 
     -        Cval(cc12,iint2(Gen4)))*Mf2(3,Gen4))/SW - 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*
     -        (2*EL*Cval(cc2,iint2(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*
     -        (2*EL*Cval(cc1,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (3*Alfa)/(2.D0*MW*Pi*SW2)*
     -        (EL*Cval(cc12,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (9*Alfa)/(16.D0*MW*Pi*SW2)*
     -        (2/3.D0*EL*Cval(cc2,iint4(Gen4))*Mf2(4,Gen4))/SW 

         loopa3 = loopa3 + 1/(16.D0*MW*Pi*SW2)*(Cval(cc0,iint2(Gen4))*
     -        (-(3*Alfa*(-2*EL*Mf2(3,Gen4)))))/SW - 
     -        (3*Alfa)/(8.D0*MW*Pi*SW2)*
     -        (-2*EL*Cval(cc1,iint2(Gen4))*Mf2(3,Gen4))/SW - 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*
     -        (-2*EL*Cval(cc2,iint2(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*
     -        (2*EL*Cval(cc1,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (9*Alfa)/(16.D0*MW*Pi*SW2)*
     -        (-2/3.D0*EL*Cval(cc2,iint4(Gen4))*Mf2(4,Gen4))/SW 
      end do

**    Lepton contributions
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            loopa1 = 1/(16.D0*MW*Pi*SW2)*(16*MW*Pi*SW2*(SW*loopa1) + 
     -           Alfa*((q12+3*q22-q32)*EL*Cval(cc1,iintL2(Gen4))*
     -           Mf2(2,Gen4)) + 3*Alfa*((q12-1/3.D0*q32+1/3.D0*q22)*EL*
     -           Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4)) + Alfa*(EL*
     -           (4*Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4) - 
     -           8*Cval(cc00,iintL2(Gen4))*Mf2(2,Gen4))))/SW
            
            loopa2 = 1/(16.D0*MW*Pi*SW2)*(16*MW*Pi*SW2*(SW*loopa2) + 
     -           Alfa*(2*EL*Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4)) + 
     -           8*Alfa*(EL*Cval(cc12,iintL2(Gen4))*Mf2(2,Gen4)) + 
     -           3*Alfa*(2/3.D0*EL*Cval(cc2,iintL2(Gen4))*
     -           Mf2(2,Gen4)))/SW
            
            loopa3 = 1/(16.D0*MW*Pi*SW2)*(16*MW*Pi*SW2*(SW*loopa3) + 
     -           Alfa*(2*EL*Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4)) + 
     -           3*Alfa*(-2/3.D0*EL*Cval(cc2,iintL2(Gen4))*
     -           Mf2(2,Gen4)))/SW
         end do
      end if


** Now we want the fermion form factors
      a2vtxWfm = -loopa2
      a1vtxWfm = loopa1 + (loopa2*q1q2) + ctFM       
      a3vtxWfm = loopa3*(0.D0,-1.D0)


** Bosonic contributions
      if (sector .eq. 3) then
         loopa1 = loopa1 + (Alfa*MW*(1 - CW2 - SW2))/(2.D0*Pi*SW2)*(EL)/ 
     -        SW + Cval(cc0,iintB14)*((5*Alfa*MW)/(8.D0*Pi*SW2)*
     -        ((q22+q12-q32)*EL)/SW + 
     -        (Alfa*MW*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2**2*Pi*SW2)*
     -        EL/SW) - (3*Alfa)/(32.D0*MW*Pi*SW2)*
     -        (EL*MH2*Bval(bb0,iintB1))/SW + Alfa/(16.D0*MW*Pi*SW2)*
     -        (EL*(-MH2 + MW2*(-12 + 12*SW2))*Bval(bb0,iintB2))/SW -
     -        Alfa/(32.D0*MW*Pi*SW2)*
     -        (EL*(MH2 + 8*MW2)*Bval(bb0,iintB3))/SW + (Alfa*MW*
     -        (CW2*(4*CW2 - 2*SW2) - SW2**2))/(8.D0*CW2*Pi*SW2)*
     -        (EL*Bval(bb0,iintB4))/SW - (Alfa*MW)/(8.D0*Pi)*
     -        (EL*(Bval(bb0,iintB5) + Bval(bb0,iintB6)))/SW - 
     -        (Alfa*MW)/(8.D0*Pi*SW2)*
     -        (EL*(Bval(bb0,iintB7) + Bval(bb0,iintB8)))/SW - 
     -        (Alfa*MW)/(8.D0*CW2*Pi)*(EL*SW*Bval(bb0,iintB9)) - 
     -        (3*Alfa*MW)/(8.D0*Pi*SW2)*(EL*MH2*Cval(cc0,iintB10))/SW - 
     -        (Alfa*MW*MW2)/(4.D0*Pi*SW2)*(EL*Cval(cc0,iintB13))/SW + 
     -        (3*Alfa)/(8.D0*MW*Pi*SW2)*
     -        (EL*MH2*Cval(cc00,iintB10))/SW + Alfa/(8.D0*MW*Pi*SW2)*
     -        (EL*(MH2 + 18*CW2*MW2 +2*MW2*SW2)*Cval(cc00,iintB11))/SW + 
     -        (2*Alfa*MW)/Pi*(EL*Cval(cc00,iintB12))/SW + 
     -        Alfa/(8.D0*MW*Pi*SW2)*
     -        (EL*(MH2 + 2*MW2)*Cval(cc00,iintB13))/SW + 
     -        Cval(cc0,iintB11)*((Alfa*MW*(2*CW2 - SW2))/(8.D0*Pi*SW2)*
     -        (q22*EL)/SW + (Alfa*MW*(10*CW2 + SW2))/(16.D0*Pi*SW2)*
     -        ((q22+q12-q32)*EL)/SW - (Alfa*MW*(2*CW2 + SW2))/
     -        (16.D0*Pi*SW2)*((q22+q12-q32)*EL)/SW + (Alfa*MW)/
     -        (8.D0*Pi)*((q12-2*q32)*EL)/SW+(Alfa*MW)/(8.D0*CW2*Pi*SW2)*
     -        (EL*(CW2*(4*CW2*MW2 - 2*MW2*SW2) - SW2**2*MH2))/
     -        SW) + Alfa/(8.D0*CW2*MW*Pi*SW2)*(EL*Cval(cc00,iintB14)*
     -        (-2*CW*CW2*MW*MZ + 2*MW2*SW2 + CW2*(MH2 + 20*MW2)))/SW +
     -        1/(8.D0*Pi)*(Cval(cc0,iintB12)*(4*Alfa*MW*
     -        ((q12+q22-3/4.D0*q32)*EL) - Alfa*MW*(EL*MH2)))/SW

         loopa1 = loopa1 + Cval(cc2,iintB11)*((Alfa*MW*(5*CW2 - 4*SW2))/
     -        (8.D0*Pi*SW2)*(q22*EL)/SW - 
     -        (Alfa*MW*(2*CW2 + SW2))/(16.D0*Pi*SW2)*((q12-q32)*EL)/SW + 
     -        (Alfa*MW)/(16.D0*Pi)*((q12-q32)*EL)/SW) + 
     -        Cval(cc1,iintB11)*(-((Alfa*MW)/(8.D0*Pi)*(q32*EL)/SW) + 
     -        (3*Alfa*MW*(-2*CW2 + SW2))/(16.D0*Pi*SW2)*((q12-q22-q32)*
     -        EL)/SW + (Alfa*MW*(2*CW2 + SW2))/(16.D0*Pi*SW2)*
     -        ((q32+q12-q22)*EL)/SW) + (Alfa*MW)/(8.D0*Pi)*
     -        ((3*q22+q12-q32)*EL*Cval(cc1,iintB12))/SW + 
     -        (3*Alfa*MW)/(8.D0*Pi)*((q12-1/3.D0*q32+1/3.D0*q22)*EL*
     -        Cval(cc2,iintB12))/SW + 1/(8.D0*CW2*Pi*SW2)*
     -        (Cval(cc1,iintB14)*(Alfa*CW2*MW*((3*q22+q12-q32)*EL)))/SW+ 
     -        1/(8.D0*CW2*Pi*SW2)*(Cval(cc2,iintB14)*
     -        (3*Alfa*CW2*MW*((q12-1/3.D0*q32+1/3.D0*q22)*EL)))/SW


         loopa2 = loopa2 + Cval(cc0,iintB14)*((5*Alfa*MW)/(8.D0*Pi*SW2)*
     -        (8/5.D0*EL)/SW) + Cval(cc0,iintB11)*((Alfa*MW*
     -        (CW2 - SW2))/(4.D0*Pi*SW2)*(EL)/SW - 
     -        (Alfa*MW*(2*CW2 + SW2))/(16.D0*Pi*SW2)*
     -        (-8*EL)/SW) + Cval(cc2,iintB11)*((Alfa*MW*(CW2 - SW2))/
     -        (4.D0*Pi*SW2)*EL/SW + 
     -        (Alfa*MW)/(16.D0*Pi)*(-8*EL)/SW) + 
     -        Cval(cc1,iintB11)*(Alfa/(8.D0*MW*Pi*SW2)*
     -        (EL*(MH2 + 18*CW2*MW2 + 2*MW2*SW2))/SW) + 
     -        (Alfa*MW)/(8.D0*Pi)*(-8*EL*Cval(cc1,iintB12))/SW + 
     -        (3*Alfa)/(8.D0*MW*Pi*SW2)*
     -        (EL*MH2*(Cval(cc1,iintB10) + Cval(cc11,iintB10) + 
     -        Cval(cc12,iintB10)))/SW + Alfa/(8.D0*MW*Pi*SW2)*
     -        (EL*(MH2 + 18*CW2*MW2 + 2*MW2*SW2)*
     -        (Cval(cc11,iintB11) + Cval(cc12,iintB11)))/SW - 
     -        (2*Alfa*MW)/Pi*(EL*Cval(cc12,iintB12))/SW - Alfa/
     -        (8.D0*MW*Pi*SW2)*(EL*(MH2 + 2*MW2)*Cval(cc12,iintB13))/SW+ 
     -        (3*Alfa*MW)/(8.D0*Pi)*(-8/3.D0*EL*Cval(cc2,iintB12))/SW - 
     -        (Alfa*MW)/(4.D0*Pi*SW2)*
     -        (EL*(Cval(cc1,iintB13) + Cval(cc2,iintB13)))/SW + 
     -        Alfa/(8.D0*CW2*MW*Pi*SW2)*(EL*Cval(cc12,iintB14)*
     -        (2*CW*CW2*MW*MZ - 2*MW2*SW2 + CW2*(-MH2 - 20*MW2)))/SW +
     -        1/(8.D0*CW2*Pi*SW2)*(Cval(cc1,iintB14)*
     -        (Alfa*((-(2*CW2*MW*EL) - 2*MW*SW2*EL))))/SW + 
     -        1/(8.D0*CW2*Pi*SW2)*(Cval(cc2,iintB14)*
     -        (Alfa*((-(2*CW2*MW*EL) - 2*MW*SW2*EL))))/SW         
      end if


** This converts between FC notation and standard notation (in FC 
** Eps = -i*levi-civita) for the a3 part of the loop
      loopa3 = loopa3*(0.D0,-1.D0)


** Debugging output
      if (debug .ne. 0) then
         write(*,*)'WWH loopa1 =', loopa1
         write(*,*)'WWH loopa2 =', loopa2
         write(*,*)'WWH loopa3 =', loopa3
      end if


** Now we want the form factors themselves 
      a2vtxW = a2vtxW - loopa2
      
      a1vtxW = a1vtxW + loopa1 + (loopa2*q1q2) + ct 
      
      a3vtxW = a3vtxW + loopa3

      if (debug .ne. 0) then
         write(*,*)"WWH vtx a1 =", a1vtxW
         write(*,*)"WWH vtx a2 =", a2vtxW
         write(*,*)"WWH vtx a3 =", a3vtxW
         write(*,*)"WWH FERMION contribution to vtx a1 =", a1vtxWfm
         write(*,*)"WWH FERMION contribution to vtx a2 =", a2vtxWfm
         write(*,*)"WWH FERMION contribution to vtx a3 =", a3vtxWfm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the ZZ-higgs vertex
*** in the standard model

      subroutine zz_vtxSM(a1vtxZ,a2vtxZ,a3vtxZ,q12,q22,q32,q1q2,debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** The arguments of the integrals
      integer*8 iint1(3), iint2(3), iint3(3), iint4(3)
      integer*8 iintL1(3), iintL2(3)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6, iintB7, iintB8
      integer*8 iintB9, iintB10

** These are the contributions from the counter-term (g-mu-nu), the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1, loopa2, loopa3


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for ZZ-higgs 
      double complex a1vtxZ, a2vtxZ, a3vtxZ


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


** Dummy parameter
      double precision dummy


** Particle generation do loop variables
      integer g4ini, Gen4


      double precision Alfa, Alfa2

      integer debug


** This is needed for looptools 
#include "looptools.h"

** Declaration of parameters
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions and 3 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


** Now for the arguments of the integrals - first the quarks arguments
      do Gen4 = g4ini, 3         
         iint1(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))         
         iint2(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),
     -        Mf2(3,Gen4),Mf2(3,Gen4))         
         iint3(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))         
         iint4(Gen4) = Cget(q22,q32,q12,Mf2(4,Gen4),
     -        Mf2(4,Gen4),Mf2(4,Gen4))
      enddo 

** Now the leptons
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))         
            iintL2(Gen4) = Cget(q22,q32,q12,Mf2(2,Gen4),
     -           Mf2(2,Gen4),Mf2(2,Gen4))
         enddo
      end if

** Now the bosons
      if (sector .eq. 3) then
         iintB1 = Bget(q12,MH2,MZ2)
         iintB2 = Bget(q22,MH2,MZ2)
         iintB3 = Bget(q12,MW2,MW2)
         iintB4 = Bget(q22,MW2,MW2)
         iintB5 = Bget(q32,MH2,MH2)
         iintB6 = Bget(q32,MW2,MW2)
         iintB7 = Bget(q32,MZ2,MZ2)
         iintB8 = Cget(q22,q32,q12,MH2,MZ2,MZ2)
         iintB9 = Cget(q22,q32,q12,MW2,MW2,MW2)
         iintB10 = Cget(q32,q12,q22,MH2,MH2,MZ2)
      end if


*** First, reset each component to zero, just in case
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0) 
      ct = (0D0,0D0)   


** Now the insanely long and hideously complicated form factors.  In bits.

** First the counterterm 
      ct = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(4*MW2*SW2*(dSW1*SW) + 
     -     CW2*(SW2*dMWsq1 + 2*MW2*SW2*dZet1 + MW2*SW2*dZH1 + 
     -     2*MW2*SW2*dZZZ1 - 2*MW2*(dSW1*SW))))/SW

** counterterm for the fermions
      ctFM = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(4*MW2*SW2*(dSW1fm*SW) + 
     -     CW2*(SW2*dMWsq1fm + 2*MW2*SW2*dZet1fm + MW2*SW2*dZH1fm + 
     -     2*MW2*SW2*dZZZ1fm - 2*MW2*(dSW1fm*SW))))/SW
            
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'ZZH ct=', ct
         write(*,*)'ZZH fermion ct=', ctFM
      end if


** Now the loop
 
** For the quarks
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + (Alfa*(9 + SW2*(-12 + 16*SW2)))/(12.D0*CW2*
     -        MW*Pi*SW2)*(EL*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4))/SW + 
     -        (Alfa*(-9 + (24 - 32*SW2)*SW2))/(6.D0*CW2*MW*Pi*SW2)*
     -        (EL*Cval(cc00,iint2(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -        ((q12+3*q22-q32)*EL*Cval(cc1,iint2(Gen4))*Mf2(3,Gen4))/SW+ 
     -        (9*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -        ((q12-1/3.D0*q32+1/3.D0*q22)*EL*Cval(cc2,iint2(Gen4))*
     -        Mf2(3,Gen4))/SW + Cval(cc0,iint2(Gen4))*
     -        ((Alfa*(1/2.D0 - (2*SW2)/3.D0))/(CW2*MW*Pi)*
     -        ((q12+q22-q32)*EL*Mf2(3,Gen4))/SW + 
     -        (3*Alfa)/(4.D0*CW2*MW*Pi*SW2)*(EL*Mf2(3,Gen4)**2)/SW) + 
     -        (Alfa*(9 + SW2*(-6 + 4*SW2)))/(12.D0*CW2*MW*Pi*SW2)*
     -        (EL*Bval(bb0,iint3(Gen4))*Mf2(4,Gen4))/SW + 
     -        (Alfa*(-9 + (12 - 8*SW2)*SW2))/(6.D0*CW2*MW*Pi*SW2)*
     -        (EL*Cval(cc00,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (3*Alfa)/(16.D0*CW2*MW*Pi*SW2)*((q12+3*q22-q32)*
     -        EL*Cval(cc1,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (9*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -        ((q12-1/3.D0*q32+1/3.D0*q22)*EL*Cval(cc2,iint4(Gen4))*
     -        Mf2(4,Gen4))/SW + Cval(cc0,iint4(Gen4))*
     -        ((Alfa*(1/2.D0 - SW2/3.D0))/(2.D0*CW2*MW*Pi)*
     -        ((q12+q22-q32)*EL*Mf2(4,Gen4))/SW + 
     -        (3*Alfa)/(4.D0*CW2*MW*Pi*SW2)*(EL*Mf2(4,Gen4)**2)/SW)

         loopa2 = loopa2 + (3*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -        (2*EL*Cval(cc1,iint2(Gen4))*Mf2(3,Gen4))/SW + 
     -        (Alfa*(9 + SW2*(-24 + 32*SW2)))/(6.D0*CW2*MW*Pi*SW2)*
     -        (EL*Cval(cc12,iint2(Gen4))*Mf2(3,Gen4))/SW + 
     -        (9*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -        (2/3.D0*EL*Cval(cc2,iint2(Gen4))*Mf2(3,Gen4))/SW + 
     -        Cval(cc0,iint2(Gen4))*((Alfa*(1/2.D0 - (2*SW2)/
     -        3.D0))/(CW2*MW*Pi)*(2*EL*Mf2(3,Gen4))/SW) + 
     -        (3*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -        (2*EL*Cval(cc1,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (Alfa*(9 + SW2*(-12 + 8*SW2)))/(6.D0*CW2*MW*Pi*SW2)*
     -        (EL*Cval(cc12,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        (9*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -        (2/3.D0*EL*Cval(cc2,iint4(Gen4))*Mf2(4,Gen4))/SW + 
     -        Cval(cc0,iint4(Gen4))*((Alfa*(1/2.D0 - SW2/3.D0))/
     -        (2.D0*CW2*MW*Pi)*(2*EL*Mf2(4,Gen4))/SW)
      end do


** For the leptons
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            loopa1 = loopa1 + 
     -           (Alfa*(1 + SW2*(-2 + 4*SW2)))/(4.D0*CW2*MW*Pi*SW2)*
     -           (EL*Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4))/SW + 
     -           (Alfa*(-1 + (4 - 8*SW2)*SW2))/(2.D0*CW2*MW*Pi*SW2)*
     -           (EL*Cval(cc00,iintL2(Gen4))*Mf2(2,Gen4))/SW + 
     -           Alfa/(16.D0*CW2*MW*Pi*SW2)*((q12+3*q22-q32)*
     -           EL*Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/SW + 
     -           (3*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -           ((q12-1/3.D0*q32+1/3.D0*q22)*EL*Cval(cc2,iintL2(Gen4))*
     -           Mf2(2,Gen4))/SW + Cval(cc0,iintL2(Gen4))*
     -           ((Alfa*(1/2.D0 - SW2))/(2.D0*CW2*MW*Pi)*
     -           ((q12+q22-q32)*EL*Mf2(2,Gen4))/SW + 
     -           Alfa/(4.D0*CW2*MW*Pi*SW2)*(EL*Mf2(2,Gen4)**2)/SW)

            loopa2 = loopa2 + Alfa/(16.D0*CW2*MW*Pi*SW2)*
     -           (2*EL*Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/SW + 
     -           (Alfa*(1 + SW2*(-4 + 8*SW2)))/(2.D0*CW2*MW*Pi*SW2)*
     -           (EL*Cval(cc12,iintL2(Gen4))*Mf2(2,Gen4))/SW + 
     -           (3*Alfa)/(16.D0*CW2*MW*Pi*SW2)*
     -           (2/3.D0*EL*Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4))/SW + 
     -           Cval(cc0,iintL2(Gen4))*((Alfa*(1/2.D0 - SW2))/
     -           (2.D0*CW2*MW*Pi)*(2*EL*Mf2(2,Gen4))/SW)
         end do
      end if

** form factors for the fermions 
      a2vtxZfm = -loopa2
      a1vtxZfm = loopa1 + (loopa2*q1q2) + ctFM
      a3vtxZfm = (0D0,0D0)


**    The bosonic contributions  
      if (sector .eq. 3) then
         loopa1 = loopa1 + Cval(cc1,iintB9)*((3*Alfa*CW2*MW)/
     -        (4.D0*Pi*SW2)*(q22*EL)/SW + 
     -        (Alfa*CW2*MW)/(4.D0*Pi*SW2)*((q12-q32)*EL)/SW) - 
     -        (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -        (EL*(Bval(bb0,iintB1) + Bval(bb0,iintB2)))/SW - 
     -        (Alfa*MW)/(4.D0*CW2*Pi)*
     -        (EL*SW*(Bval(bb0,iintB3) + Bval(bb0,iintB4))) - 
     -        (3*Alfa)/(32.D0*CW2*MW*Pi*SW2)*(EL*MH2*Bval(bb0,iintB5))/
     -        SW - Alfa/(32.D0*CW2*MW*Pi*SW2)*(EL*MH2*Bval(bb0,iintB7))/
     -        SW - (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(EL*MH2*
     -        Cval(cc0,iintB10))/SW -(Alfa*MW*MW2)/(4.D0*CW2**3*Pi*SW2)*
     -        (EL*Cval(cc0,iintB8))/SW + (3*Alfa)/(8.D0*CW2*MW*Pi*SW2)*
     -        (EL*MH2*Cval(cc00,iintB10))/SW + Alfa/(8.D0*CW2**2*
     -        MW*Pi*SW2)*(EL*Cval(cc00,iintB8)*(2*MW2 + CW2*MH2))/SW + 
     -        Cval(cc0,iintB9)*((Alfa*MW*(-10*CW2 + SW2))/(8.D0*Pi*SW2)*
     -        (q32*EL)/SW + (Alfa*MW*(5*CW2 - SW2))/(4.D0*Pi*SW2)*
     -        ((q12+q22)*EL)/SW + (3*Alfa*MW)/(8.D0*Pi)*(q32*EL)/SW + 
     -        (Alfa*MW)/(4.D0*CW2*Pi*SW2)*(EL*(-2*MW2*SW2**2 + 
     -        CW2*(4*CW2*MW2 + 2*MW2*SW2) - SW2**2*MH2))/SW) + Alfa/
     -        (16.D0*CW2*MW*Pi*SW2)*(EL*Bval(bb0,iintB6)*
     -        (-(SW2**2*MH2) + CW2*(8*MW2*SW2 + 2*SW2*MH2 + CW2*
     -        (-MH2 - 8*MW2))))/SW + Alfa/(4.D0*CW2*MW*Pi*SW2)*
     -        (EL*Cval(cc00,iintB9)*(2*MW2*SW2**2 + SW2**2*MH2 + 
     -        CW2*(-4*MW2*SW2 - 2*SW2*MH2 + CW2*(MH2 + 18*MW2))))/SW + 
     -        Cval(cc2,iintB9)*((3*Alfa*CW2*MW)/(4.D0*Pi*SW2)*(q12*EL)/ 
     -        SW + (Alfa*CW2*MW)/(4.D0*Pi*SW2)*((q22-q32)*EL)/SW)


         loopa2 = loopa2 + Cval(cc1,iintB9)*(-(Alfa*MW*
     -        (SW2**2 + CW2*(CW2 + 2*SW2)))/(2.D0*CW2*Pi*SW2)*EL/SW) + 
     -        (3*Alfa)/(8.D0*CW2*MW*Pi*SW2)*
     -        (EL*MH2*(Cval(cc1,iintB10) + Cval(cc11,iintB10) + 
     -        Cval(cc12,iintB10)))/SW - 
     -        (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -        (EL*(Cval(cc1,iintB8) + Cval(cc2,iintB8)))/SW - 
     -        Alfa/(8.D0*CW2**2*MW*Pi*SW2)*
     -        (EL*Cval(cc12,iintB8)*(2*MW2 + CW2*MH2))/SW + 
     -        Cval(cc0,iintB9)*((Alfa*MW*(4*CW2 - SW2))/(2.D0*Pi*SW2)* 
     -        EL/SW + (3*Alfa*MW)/(8.D0*Pi)*(-4*EL)/SW) + 
     -        Alfa/(4.D0*CW2*MW*Pi*SW2)*(EL*Cval(cc12,iintB9)*
     -        (-2*MW2*SW2**2 - SW2**2*MH2 + 
     -        CW2*(4*MW2*SW2 + 2*SW2*MH2 + CW2*(-MH2 - 18*MW2))))/
     -        SW + Cval(cc2,iintB9)*(-(Alfa*MW*(SW2**2 + CW2*
     -        (CW2 + 2*SW2)))/(2.D0*CW2*Pi*SW2)*EL/SW)
      end if


** The eps part of the loop 
      loopa3 = (0.D0,0.D0)

      if (debug .ne. 0) then
         write(*,*)'ZZH loopa1=', loopa1
         write(*,*)'ZZH loopa2=', loopa2
         write(*,*)'ZZH loopa3 =', loopa3
      end if


** Now we want the form factors themselves 
      a2vtxZ = a2vtxZ - loopa2

      a1vtxZ = a1vtxZ + loopa1 + (loopa2*q1q2) + ct 

      a3vtxZ = a3vtxZ + loopa3

      if (debug .ne. 0) then
         write(*,*)"ZZH vtx a1 =", a1vtxZ
         write(*,*)"ZZH vtx a2 =", a2vtxZ
         write(*,*)"ZZH vtx a3 =", a3vtxZ
         write(*,*)"ZZH FERMION contribution to vtx a1 =", a1vtxZfm
         write(*,*)"ZZH FERMION contribution to vtx a2 =", a2vtxZfm
         write(*,*)"ZZH FERMION contribution to vtx a3 =", a3vtxZfm
      end if
 

      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the photon-photon-H 
*** vertex

      subroutine aa_vtxSM(a1vtxAH,a2vtxAH,a3vtxAH,q12,q22,q32,q1q2,
     -     debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** The arguments of the integrals
      integer*8 iint1(3), iint2(3), iint3(3), iint4(3)
      integer*8 iintL1(3), iintL2(3)
      integer*8 iintB1, iintB2, iintB3, iintB4

** These are the contributions from the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex loopa1, loopa2, loopa3


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for photon-photon-higgs 
      double complex a1vtxAH, a2vtxAH, a3vtxAH


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions
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


** Particle do loop parameters
      integer g4ini, Gen4

      double precision Alfa, Alfa2

      integer debug


** This is needed for looptools. 
#include "looptools.h"

** Declaration of parameters
#include "mssm.inc"
#include "global.inc"


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions and 3 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Now for the arguments of the integrals - first the fermionic arguments, 
** then the bosonic. 
      do Gen4 = g4ini, 3   
         iint1(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
         iint2(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),Mf2(3,Gen4),
     -        Mf2(3,Gen4))
         iint3(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
         iint4(Gen4) = Cget(q22,q32,q12,Mf2(4,Gen4),Mf2(4,Gen4),
     -        Mf2(4,Gen4))   
      enddo 

      if (sector .ne. 1) then
         do Gen4 = 1, 3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL2(Gen4) = Cget(q22,q32,q12,Mf2(2,Gen4),Mf2(2,Gen4),
     -           Mf2(2,Gen4))
         enddo
      end if

      if (sector .eq. 3) then
         iintB1 = Bget(q32,MW2,MW2)
         iintB2 = Bget(q22,MW2,MW2)
         iintB3 = Bget(q12,MW2,MW2)
         iintB4 = Cget(q22,q32,q12,MW2,MW2,MW2)
      end if

*** First, reset to zero,just in case
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0) 


** Now the insanely long and hideously complicated form factors.  In bits.

** Now the loop
** For the quarks
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + 1/(6.D0*MW*Pi)*(Alfa*((q12+q22-q32)*EL*(-4*
     -        Cval(cc0,iint2(Gen4))*Mf2(3,Gen4) - Cval(cc0,iint4(Gen4))*
     -        Mf2(4,Gen4))) + Alfa*(EL*(8*Bval(bb0,iint1(Gen4))*
     -        Mf2(3,Gen4) - 32*Cval(cc00,iint2(Gen4))*Mf2(3,Gen4) + 
     -        2*Bval(bb0,iint3(Gen4))*Mf2(4,Gen4) - 8*
     -        Cval(cc00,iint4(Gen4))*Mf2(4,Gen4))))/SW

         loopa2 = loopa2 + 1/(6.D0*MW*Pi)*(Alfa*(2*EL*(-4*
     -        Cval(cc0,iint2(Gen4))*Mf2(3,Gen4) - Cval(cc0,iint4(Gen4))*
     -        Mf2(4,Gen4))) + Alfa*(EL*(32*Cval(cc12,iint2(Gen4))*
     -        Mf2(3,Gen4) + 8*Cval(cc12,iint4(Gen4))*Mf2(4,Gen4))))/SW
      end do

** For the leptons
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            loopa1 = loopa1 + 1/(2.D0*MW*Pi)*(-(Alfa* ((q12+q22-q32)*EL*
     -           Cval(cc0,iintL2(Gen4))*Mf2(2,Gen4))) + Alfa*(EL*(2*
     -           Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4) - 
     -           8*Cval(cc00,iintL2(Gen4))*Mf2(2,Gen4))))/SW

            loopa2 = loopa2 + 1/(2.D0*MW*Pi)*(-(Alfa*(2*EL*
     -           Cval(cc0,iintL2(Gen4))*Mf2(2,Gen4))) + 8*Alfa*(EL*
     -           Cval(cc12,iintL2(Gen4))*Mf2(2,Gen4)))/SW
         end do
      end if


** form factors for the fermions
      a2vtxAHfm = -loopa2
      a1vtxAHfm = loopa1 + (loopa2*q1q2)  
      a3vtxAHfm = (0D0,0D0)


** The bosonic contributions  
      if (sector .eq. 3) then
         loopa1 = loopa1 - Alfa/(4.D0*MW*Pi)*(EL*(MH2 + 4*MW2)*
     -        Bval(bb0,iintB1))/SW - (Alfa*MW)/(4.D0*Pi)*(EL*
     -        (Bval(bb0,iintB2) + Bval(bb0,iintB3)))/SW + Alfa/(MW*Pi)*
     -        (EL*(MH2 + 6*MW2)*Cval(cc00,iintB4))/SW + (Alfa*MW)/(4.D0*
     -        Pi)*((3*q22-q32+q12)*EL*Cval(cc1,iintB4))/SW + (Alfa*MW)/
     -        (4.D0*Pi)*((3*(q12+1/3.D0*q22-1/3.D0*q32) -(q12-q22+q32) + 
     -        (q12-q22+q32))*EL*Cval(cc2,iintB4))/SW + 1/(4.D0*Pi)*
     -        (Cval(cc0,iintB4)*(5*Alfa*MW*((q12+q22-q32)*EL) + Alfa*MW*
     -        ((q12-q32)*EL) + Alfa*MW*((q22-q32)*EL) - Alfa*MW*
     -        (EL*MH2)))/SW

         loopa2 = loopa2 - Alfa/(MW*Pi)*(EL*(MH2 + 6*MW2)*
     -        Cval(cc12,iintB4))/SW + 1/(4.D0*Pi)*(Cval(cc0,iintB4)*
     -        (5*Alfa*MW*((8/5.D0)*EL) + Alfa*MW*4*EL +Alfa*MW*4*EL))/SW
      end if


** The eps part of the loop 
      loopa3 = (0.D0,0.D0)

      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'AAH loopa1=', loopa1
         write(*,*)'AAH loopa2=', loopa2
         write(*,*)'AAH loopa3 =', loopa3
      end if


** Now we want the form factors themselves 
      a2vtxAH = 0D0 - loopa2

      a1vtxAH = loopa1 + (loopa2*q1q2)  

      a3vtxAH = loopa3

      if (debug .ne. 0) then
         write(*,*)"AAH vtx a1 =", a1vtxAH
         write(*,*)"AAH vtx a2 =", a2vtxAH
         write(*,*)"AAH vtx a3 =", a3vtxAH
         write(*,*)"AAH FERMION contribution to vtx a1 =", a1vtxAHfm
         write(*,*)"AAH FERMION contribution to vtx a2 =", a2vtxAHfm
         write(*,*)"AAH FERMION contribution to vtx a3 =", a3vtxAHfm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the AZH vertex

      subroutine az_vtxSM(a1vtxAZH,a2vtxAZH,a3vtxAZH,q12,q22,q32,q1q2,
     -     debug,fm)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** The arguments of the integrals
      integer*8 iint1(3), iint2(3), iint3(3), iint4(3)
      integer*8 iintL1(3), iintL2(3)
      integer*8 iintB1, iintB2, iintB3, iintB4

** These are the contributions from the counter-term (g-mu-nu), the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1, loopa2, loopa3

** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for photon-Z-higgs 
      double complex a1vtxAZH, a2vtxAZH, a3vtxAZH


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions.  'fm' = 1 for AZ and 2 for AZ2
      integer fm
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


** Particle do loop parameters
      integer g4ini, Gen4

      
      double precision Alfa, Alfa2

      integer debug


** This is needed for looptools. 
#include "looptools.h"

** Declaration of parameters
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions and 3 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


** Now for the arguments of the integrals - first the fermionic arguments, 
** then the bosonic. 
      do Gen4 = g4ini, 3  
         iint1(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
         iint2(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),Mf2(3,Gen4),
     -        Mf2(3,Gen4))
         iint3(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
         iint4(Gen4) = Cget(q22,q32,q12,Mf2(4,Gen4),Mf2(4,Gen4),
     -        Mf2(4,Gen4))
      enddo 

      if (sector .ne. 1) then
         do Gen4 = 1, 3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL2(Gen4) = Cget(q22,q32,q12,Mf2(2,Gen4),Mf2(2,Gen4),
     -           Mf2(2,Gen4))
         enddo
      end if

      if (sector .eq. 3) then
         iintB1 = Bget(q32,MW2,MW2)
         iintB2 = Bget(q22,MW2,MW2)
         iintB3 = Bget(q12,MW2,MW2)
         iintB4 = Cget(q22,q32,q12,MW2,MW2,MW2)
      end if


*** First, reset to zero,just in case
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0) 
      ct = (0D0,0D0)   


** Now the insanely long and hideously complicated form factors.  In bits.

** First the counterterm - note this is only calculated on the final run 
** through of the fermions
      ct = MW/(2.D0*CW2)*(dZZA1*EL)/SW

** counterterm for fermions
      ctFM = MW/(2.D0*CW2)*(dZZA1fm*EL)/SW
      
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'AZH ct=', ct
         write(*,*)'AZH fermion ct=', ctFM
      end if


** Now the loop
      
** For the quarks
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + 1/(24.D0*CW*MW*Pi*SW2)*(Alfa*((q12+q22-q32)*
     -        EL*(SW2*(-16*Cval(cc0,iint2(Gen4))*Mf2(3,Gen4) - 4*
     -        Cval(cc0,iint4(Gen4))*Mf2(4,Gen4)) + 6*
     -        Cval(cc0,iint2(Gen4))*Mf2(3,Gen4) + 3*
     -        Cval(cc0,iint4(Gen4))*Mf2(4,Gen4))) + Alfa*(EL*(SW2*(32*
     -        Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) - 
     -        128*Cval(cc00,iint2(Gen4))*Mf2(3,Gen4) + 
     -        8*Bval(bb0,iint3(Gen4))*Mf2(4,Gen4) - 
     -        32*Cval(cc00,iint4(Gen4))*Mf2(4,Gen4)) - 
     -        12*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) + 
     -        48*Cval(cc00,iint2(Gen4))*Mf2(3,Gen4) - 
     -        6*Bval(bb0,iint3(Gen4))*Mf2(4,Gen4) + 
     -        24*Cval(cc00,iint4(Gen4))*Mf2(4,Gen4))))

         loopa2 = loopa2 + 1/(24.D0*CW*MW*Pi*SW2)*(Alfa*(2*EL*
     -        (SW2*(-16*Cval(cc0,iint2(Gen4))*Mf2(3,Gen4) - 
     -        4*Cval(cc0,iint4(Gen4))*Mf2(4,Gen4)) + 
     -        6*Cval(cc0,iint2(Gen4))*Mf2(3,Gen4) + 
     -        3*Cval(cc0,iint4(Gen4))*Mf2(4,Gen4))) + Alfa*(EL*(SW2*
     -        (128*Cval(cc12,iint2(Gen4))*Mf2(3,Gen4) + 
     -        32*Cval(cc12,iint4(Gen4))*Mf2(4,Gen4)) - 
     -        48*Cval(cc12,iint2(Gen4))*Mf2(3,Gen4) - 
     -        24*Cval(cc12,iint4(Gen4))*Mf2(4,Gen4))))
      end do

** For the leptons
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            loopa1 = loopa1+ 1/(8.D0*CW*MW*Pi*SW2)*(Alfa*((q12-q22-q32)*
     -           EL*(-(4*SW2*(Cval(cc0,iintL2(Gen4))*Mf2(2,Gen4))) + 
     -           Cval(cc0,iintL2(Gen4))*Mf2(2,Gen4))) + Alfa*(EL*
     -           (SW2*(8*Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4) - 
     -           32*Cval(cc00,iintL2(Gen4))*Mf2(2,Gen4)) - 
     -           2*Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4) + 
     -           8*Cval(cc00,iintL2(Gen4))*Mf2(2,Gen4))) )

            loopa2 = loopa2 + 1/(8.D0*CW*MW*Pi*SW2)*(Alfa*(2*EL*
     -           (-(4*SW2*(Cval(cc0,iintL2(Gen4))*Mf2(2,Gen4))) + 
     -           Cval(cc0,iintL2(Gen4))*Mf2(2,Gen4))) + Alfa*(EL*(32*
     -           SW2*(Cval(cc12,iintL2(Gen4))*Mf2(2,Gen4)) - 
     -           8*Cval(cc12,iintL2(Gen4))*Mf2(2,Gen4))))
         end do
      end if


** form factors for the fermions
      if (fm .eq. 1) then
         a2vtxAZHfm = -loopa2
         a1vtxAZHfm = loopa1 + (loopa2*q1q2) + ctFM 
         a3vtxAZHfm = (0D0,0D0)
      else
         a2vtxAZH2fm = -loopa2
         a1vtxAZH2fm = loopa1 + (loopa2*q1q2) + ctFM
         a3vtxAZH2fm = (0D0,0D0)
      end if


** The bosonic contributions  
      if (sector .eq. 3) then
         loopa1 = loopa1 + Bval(bb0,iintB1)*(-((Alfa*MW*(-3*CW2 + SW2))/
     -        (4.D0*CW*Pi*SW2)*EL) - (Alfa*(-CW2 + SW2))/(8.D0*CW*MW*Pi*
     -        SW2)*(EL*MH2))+Cval(cc00,iintB4)*((Alfa*MW*(-19*CW2 + 3*
     -        SW2))/(4.D0*CW*Pi*SW2)*EL + (Alfa*(-CW2 + SW2))/(4.D0*CW*
     -        MW*Pi*SW2)*(EL*(2*MH2 + MW2))) - (Alfa*MW)/(4.D0*CW*Pi)*
     -        (EL*(Bval(bb0,iintB2) + Bval(bb0,iintB3))) + 1/(4.D0*CW*
     -        Pi*SW2)*(Cval(cc1,iintB4)*(-Alfa*CW2*MW*((q12+3*q22-q32)*
     -        EL) + Alfa*CW2*MW*((q12-q22-q32)*EL) + Alfa*MW*SW2*((q12-
     -        q22-q32)*EL))) + Cval(cc0,iintB4)*(-((Alfa*CW*MW)/(4.D0*
     -        Pi*SW2)*((5*(q12+q22-q32) + (q12-q32))*EL)) + (Alfa*MW)/
     -        (4.D0*CW*Pi)*((q22-q32)*EL) - (Alfa*MW)/(4.D0*CW*Pi*SW2)*
     -        (EL*(3*CW2*MW2 + 3*MW2*SW2 + SW2*MH2))) + 1/(4.D0*CW*Pi*
     -        SW2)*(Cval(cc2,iintB4)*(-(3*Alfa*CW2*MW*((q12+1/3.D0*q22-
     -        1/3.D0*q32)*EL)) + Alfa*CW2*MW*((q12-q22+q32)*EL) + Alfa*
     -        MW*SW2*((q12-q22+q32)*EL)))

         loopa2 = loopa2 + Cval(cc12,iintB4)*(-((Alfa*MW*(-19*CW2 + 3*
     -        SW2))/(4.D0*CW*Pi*SW2)*EL) - (Alfa*(-CW2 + SW2))/(4.D0*CW*
     -        MW*Pi*SW2)*(EL*(2*MH2 + MW2))) + 1/(4.D0*CW*Pi*SW2)*
     -        (Cval(cc1,iintB4)*(-(4*Alfa*MW*SW2*EL) - Alfa*CW2*MW*
     -        (-2)*EL + Alfa*CW2*MW*(-2)*EL + Alfa*MW*SW2*4*EL)) + 
     -        Cval(cc0,iintB4)*(-((Alfa*CW*MW)/(4.D0*Pi*SW2)*(12*EL)) + 
     -        (Alfa*MW)/(4.D0*CW*Pi)*(4*EL)) + 1/(4.D0*CW*Pi*SW2)*
     -        (Cval(cc2,iintB4)*(-(3*Alfa*CW2*MW*(-2/3.D0)*EL) + Alfa*
     -        CW2*MW*(-4)*EL + Alfa*MW*SW2*2*EL + Alfa*((2*CW2*MW*EL - 
     -        2*MW*SW2*EL))))
      end if

** The eps part of the loop 
      loopa3 = (0.D0,0.D0)

      if (debug .ne. 0) then
         write(*,*)'AZH loopa1=', loopa1
         write(*,*)'AZH loopa2=', loopa2
         write(*,*)'AZH loopa3 =', loopa3
      end if


** Now we want the form factors themselves 
      a2vtxAZH = 0D0 - loopa2

      a1vtxAZH = loopa1 + (loopa2*q1q2) + ct 

      a3vtxAZH = loopa3

      if (debug .ne. 0) then
         write(*,*)"AZH vtx a1 =", a1vtxAZH
         write(*,*)"AZH vtx a2 =", a2vtxAZH
         write(*,*)"AZH vtx a3 =", a3vtxAZH
         write(*,*)"AZH FERMION cont to vtx a1 =", a1vtxAZHfm, a1vtxAZH2fm
         write(*,*)"AZH FERMION cont to vtx a2 =", a2vtxAZHfm, a2vtxAZH2fm
         write(*,*)"AZH FERMION cont to vtx a3 =", a3vtxAZHfm, a3vtxAZH2fm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the WW self energy
*** in the standard model

      subroutine ww_selfSM(a1selfW,q12,q22,q32,q1q2,debug)

      implicit none

** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Integral arguments
      double complex cint1(3), cintL1(3), cintB1, cintB2, cintB3
      integer*8 iint1(3), iint2(3), iintL1(3), iintL2(3)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6

** These are the contributions from the counter-term (g-mu-nu), and the loop
** (g-mu-nu).
      double complex ct, loopa1


** This is the contribution to the formfactor a1 (g-mu-nu) for WW self energy
      double complex a1selfW


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


** Fermion do loop parameters
      integer Gen4, g4ini

      double precision Alfa, Alfa2

      integer debug

** Parameter declarations
#include "mssm.inc"
#include "global.inc"

** This is needed for looptools. 
#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in vbfnlo.dat and is 1 for generation 3 quarks,
** 2 for all fermions and 3 for all particles
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


** Clearing the values of the loop and the counterterm
      loopa1 = (0D0,0D0)  
      ct = (0D0,0D0)   


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
      ct = (-(MW*(q22*dZW1*EL)) + MW*(dMWsq1*EL) + MW*MW2*(dZW1*EL))/
     -     ((-MW2 + q22)*SW) + (-(MW*(q12*dZW1*EL)) + MW*(dMWsq1*EL) + 
     -     MW*MW2*(dZW1*EL))/((-MW2 + q12)*SW)

** counterterm for the fermions
      ctFM = (-(MW*(q22*dZW1fm*EL)) + MW*(dMWsq1fm*EL) + MW*MW2*(dZW1fm*EL))/
     -     ((-MW2 + q22)*SW) + (-(MW*(q12*dZW1fm*EL)) + MW*(dMWsq1fm*EL) + 
     -     MW*MW2*(dZW1fm*EL))/((-MW2 + q12)*SW)
      
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'WW self energy ct=', ct
         write(*,*)'WW FERMION contribution to self energy ct=', ctFM
      end if


** Now the arguments of the integrals: first the quarks
      do Gen4 = g4ini, 3
         cint1(Gen4) = A0(Mf2(4,Gen4))
         iint1(Gen4) = Bget(q22,Mf2(3,Gen4),Mf2(4,Gen4))
         iint2(Gen4) = Bget(q12,Mf2(3,Gen4),Mf2(4,Gen4))
      enddo

** Now the letpons
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            cintL1(Gen4) = A0(Mf2(2,Gen4))
            iintL1(Gen4) = Bget(q22,0.D0,Mf2(2,Gen4))
            iintL2(Gen4) = Bget(q12,0.D0,Mf2(2,Gen4))
         enddo
      end if

** Now the bosons
      if (sector .eq. 3) then
         cintB1 = A0(MH2)
         cintB2 = A0(MW2)
         cintB3 = A0(MZ2)
         iintB1 = Bget(q22,0.D0,MW2)
         iintB2 = Bget(q22,MH2,MW2)
         iintB3 = Bget(q22,MW2,MZ2)
         iintB4 = Bget(q12,0.D0,MW2)
         iintB5 = Bget(q12,MH2,MW2)
         iintB6 = Bget(q12,MW2,MZ2)
      end if


** The loop component: the quarks
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + 1/(4.D0*Pi*SW2)*(3*Alfa*MW*
     -        (q22*EL*Bval(bb1,iint1(Gen4))) + Alfa*MW*(EL*
     -        (-6*Bval(bb00,iint1(Gen4)) + 3*cint1(Gen4) + 
     -        3*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4))))/
     -        ((-MW2 + q22)*SW) +  1/(4.D0*Pi*SW2)*(3*Alfa*MW*
     -        (q12*EL*Bval(bb1,iint2(Gen4))) + Alfa*MW*(EL*
     -        (-6*Bval(bb00,iint2(Gen4)) + 3*cint1(Gen4) + 
     -        3*Bval(bb0,iint2(Gen4))*Mf2(3,Gen4))))/((-MW2 + q12)*SW)
      end do

      
** The leptons
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            loopa1 = loopa1 + 1/(4.D0*Pi*SW2)*(Alfa*MW*(q22*EL*
     -           Bval(bb1,iintL1(Gen4))) + Alfa*MW*(EL*
     -           (-2*Bval(bb00,iintL1(Gen4)) + cintL1(Gen4))))/
     -           ((-MW2 + q22)*SW) + 1/(4.D0*Pi*SW2)*(Alfa*MW*
     -           (q12*EL*Bval(bb1,iintL2(Gen4))) + Alfa*MW*(EL*
     -           (-2*Bval(bb00,iintL2(Gen4)) + cintL1(Gen4))))/
     -           ((-MW2 + q12)*SW)
         end do
      end if

** formfactor for fermions
      a1selfWfm = loopa1 + ctFM

** The bosons
      if (sector .eq. 3) then
         loopa1 = loopa1 + ((Alfa*MW)/(6.D0*Pi*SW2)*(q22*EL)/SW + 
     -        (Alfa*MW*MW2*(1 - CW2 - SW2))/(2.D0*Pi*SW2)*EL/SW + 
     -        Bval(bb0,iintB3)*
     -        ((5*Alfa*CW2*MW)/(4.D0*Pi*SW2)*(q22*EL)/SW + 
     -        (Alfa*MW*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        EL/SW) - (Alfa*MW)/(16.D0*Pi*SW2)*(cintB1*EL)/SW + 
     -        (Alfa*MW*(-7 + 4*SW2))/(8.D0*Pi*SW2)*(cintB2*EL)/SW - 
     -        (Alfa*(1/4.D0 + CW2)*MW)/(4.D0*Pi*SW2)*(cintB3*EL)/SW - 
     -        (Alfa*MW*MW2)/(4.D0*Pi*SW2)*(EL*Bval(bb0,iintB2))/SW + 
     -        (2*Alfa*MW)/Pi*(EL*Bval(bb00,iintB1))/SW + 
     -        (Alfa*MW)/(4.D0*Pi*SW2)*(EL*Bval(bb00,iintB2))/SW + 
     -        (Alfa*(1/4.D0 + 2*CW2)*MW)/(Pi*SW2)*
     -        (EL*Bval(bb00,iintB3))/SW + 
     -        (Alfa*MW)/(2.D0*Pi)*(q22*EL*Bval(bb1,iintB1))/SW + 
     -        (Alfa*CW2*MW)/(2.D0*Pi*SW2)*(q22*EL*Bval(bb1,iintB3))/SW + 
     -        1/(4.D0*Pi)*(Bval(bb0,iintB1)*
     -        (5*Alfa*MW*(q22*EL) - Alfa*MW*MW2*EL))/SW)/(-MW2 + q22) + 
     -        ((Alfa*MW)/(6.D0*Pi*SW2)*(q12*EL)/SW + 
     -        (Alfa*MW*MW2*(1 - CW2 - SW2))/(2.D0*Pi*SW2)*EL/SW + 
     -        Bval(bb0,iintB6)*((5*Alfa*CW2*MW)/(4.D0*Pi*SW2)*(q12*EL)/ 
     -        SW +(Alfa*MW*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        EL/SW) - (Alfa*MW)/(16.D0*Pi*SW2)*(cintB1*EL)/SW + 
     -        (Alfa*MW*(-7 + 4*SW2))/(8.D0*Pi*SW2)*(cintB2*EL)/SW - 
     -        (Alfa*(1/4.D0 + CW2)*MW)/(4.D0*Pi*SW2)*(cintB3*EL)/SW - 
     -        (Alfa*MW*MW2)/(4.D0*Pi*SW2)*(EL*Bval(bb0,iintB5))/SW + 
     -        (2*Alfa*MW)/Pi*(EL*Bval(bb00,iintB4))/SW + 
     -        (Alfa*MW)/(4.D0*Pi*SW2)*(EL*Bval(bb00,iintB5))/SW + 
     -        (Alfa*(1/4.D0 + 2*CW2)*MW)/(Pi*SW2)*
     -        (EL*Bval(bb00,iintB6))/SW + 
     -        (Alfa*MW)/(2.D0*Pi)*(q12*EL*Bval(bb1,iintB4))/SW + 
     -        (Alfa*CW2*MW)/(2.D0*Pi*SW2)*(q12*EL*Bval(bb1,iintB6))/SW + 
     -        1/(4.D0*Pi)*(Bval(bb0,iintB4)*
     -        (5*Alfa*MW*(q12*EL) - Alfa*MW*MW2*EL))/SW)/(-MW2 + q12)
      end if
         
          
      if (debug .ne. 0) then
         write(*,*)'WW self energy loopa1=', loopa1
      end if


** Now we want the form factor itself
      a1selfW = a1selfW + loopa1 + ct


      if (debug .ne. 0) then
         write(*,*)'WW self energy a1 =', a1selfW
         write(*,*)'WW FERMION contribution to self energy a1 =', a1selfWfm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the ZZ self energy
*** in the standard model

      subroutine zz_selfSM(a1selfZ,q12,q22,q32,q1q2,debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** These are the contributions from the counter-term (g-mu-nu) and the loop
** (g-mu-nu).
      double complex ct, loopa1

** Integral arguments
      double complex cint1(3), cint2(3), cintL1(3), cintB1,cintB2,cintB3
      integer*8 iint1(3), iint2(3), iint3(3), iint4(3)
      integer*8 iintL1, iintL2, iintL3(3), iintL4(3)
      integer*8 iintB1, iintB2, iintB3, iintB4


** This is the contribution to the formfactor a1 (g-mu-nu)
      double complex a1selfZ


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


** Particle generation do loop variables
      integer g4ini, Gen4

      double precision Alfa, Alfa2

      integer debug

** This is needed for looptools. 
#include "looptools.h"

** Parameter declarations
#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions and 3 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


** Clearing the loop and counterterm values
      loopa1 = (0D0,0D0)  
      ct = (0D0,0D0)    


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
      ct = 1/CW2*(-(MW*(q22*dZZZ1*EL)) + MW*(dMZsq1*EL) + 
     -     MW*MZ2*(dZZZ1*EL))/((-MZ2 + q22)*SW) + 
     -     1/CW2*(-(MW*(q12*dZZZ1*EL)) + MW*(dMZsq1*EL) + 
     -     MW*MZ2*(dZZZ1*EL))/((-MZ2 + q12)*SW)

** counterterm for the fermions
      ctFM = 1/CW2*(-(MW*(q22*dZZZ1fm*EL)) + MW*(dMZsq1fm*EL) + 
     -     MW*MZ2*(dZZZ1fm*EL))/((-MZ2 + q22)*SW) + 
     -     1/CW2*(-(MW*(q12*dZZZ1fm*EL)) + MW*(dMZsq1fm*EL) + 
     -     MW*MZ2*(dZZZ1fm*EL))/((-MZ2 + q12)*SW)
      
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'ZZ self energy ct=', ct
         write(*,*)'ZZ FERMION contrn to self energy ct=', ctFM
      end if


** Now the loop

** Now for the arguments of the integrals: first the quarks
      do Gen4 = g4ini, 3
         cint1(Gen4) = A0(Mf2(3,Gen4))
         cint2(Gen4) = A0(Mf2(4,Gen4))
         iint1(Gen4) = Bget(q22,Mf2(3,Gen4),Mf2(3,Gen4))
         iint2(Gen4) = Bget(q22,Mf2(4,Gen4),Mf2(4,Gen4))
         iint3(Gen4) = Bget(q12,Mf2(3,Gen4),Mf2(3,Gen4))
         iint4(Gen4) = Bget(q12,Mf2(4,Gen4),Mf2(4,Gen4))
      enddo

** Now the leptons
      if (sector .ne. 1) then
         iintL1 = Bget(q22,0.D0,0.D0)
         iintL2 = Bget(q12,0.D0,0.D0)
         do Gen4 = 1, 3
            cintL1(Gen4) = A0(Mf2(2,Gen4))
            iintL3(Gen4) = Bget(q22,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL4(Gen4) = Bget(q12,Mf2(2,Gen4),Mf2(2,Gen4))
         enddo
      end if

** Now the bosons
      if (sector .eq. 3) then
         cintB1 = A0(MW2)         
         cintB2 = A0(MH2)
         cintB3 = A0(MZ2)
         iintB1 = Bget(q22,MH2,MZ2)
         iintB2 = Bget(q22,MW2,MW2)
         iintB3 = Bget(q12,MH2,MZ2)
         iintB4 = Bget(q12,MW2,MW2)
      end if

      
** The loop components: first, the quarks:
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + ((Alfa*MW*(-9 + (24 - 32*SW2)*SW2))/
     -        (12.D0*CW2**2*Pi*SW2)*(EL*Bval(bb00,iint1(Gen4)))/SW + 
     -        (Alfa*MW*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb00,iint2(Gen4)))/SW + 
     -        (Alfa*MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q22*EL*Bval(bb1,iint1(Gen4)))/SW + 
     -        (Alfa*MW*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q22*EL*Bval(bb1,iint2(Gen4)))/SW + 
     -        (Alfa*MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (EL*cint1(Gen4))/SW + 
     -        (Alfa*MW*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (EL*cint2(Gen4))/SW + (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb0,iint2(Gen4))*Mf2(4,Gen4))/SW)/(-MZ2 + q22)
         
         loopa1 = loopa1 + ((Alfa*MW*(-9 + (24 - 32*SW2)*SW2))/
     -        (12.D0*CW2**2*Pi*SW2)*(EL*Bval(bb00,iint3(Gen4)))/SW + 
     -        (Alfa*MW*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb00,iint4(Gen4)))/SW + 
     -        (Alfa*MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q12*EL*Bval(bb1,iint3(Gen4)))/SW + 
     -        (Alfa*MW*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q12*EL*Bval(bb1,iint4(Gen4)))/SW + 
     -        (Alfa*MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (EL*cint1(Gen4))/SW + 
     -        (Alfa*MW*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (EL*cint2(Gen4))/SW + (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb0,iint3(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb0,iint4(Gen4))*Mf2(4,Gen4))/SW)/(-MZ2 + q12)
      end do


** Now, the leptons ...
      if (sector .ne. 1) then
         loopa1 = loopa1 + 1/(8.D0*CW2**2*Pi*SW2)*
     -        (-(6*Alfa*MW*(EL*Bval(bb00,iintL1))) + 
     -        3*Alfa*MW*(q22*EL*Bval(bb1,iintL1)))/
     -        ((-MZ2 + q22)*SW) + 1/(8.D0*CW2**2*Pi*SW2)*
     -        (-(6*Alfa*MW*(EL*Bval(bb00,iintL2))) + 
     -        3*Alfa*MW*(q12*EL*Bval(bb1,iintL2)))/((-MZ2 + q12)*SW)

         do Gen4 = 1, 3
            loopa1 = loopa1 + ((Alfa*MW*(-1 + (4 - 8*SW2)*SW2))/
     -           (4.D0*CW2**2*Pi*SW2)*(EL*Bval(bb00,iintL3(Gen4)))/SW + 
     -           (Alfa*MW*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -           (q22*EL*Bval(bb1,iintL3(Gen4)))/SW + 
     -           (Alfa*MW*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -           (EL*cintL1(Gen4))/SW + (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -           (EL*Bval(bb0,iintL3(Gen4))*Mf2(2,Gen4))/SW)/
     -           (-MZ2 + q22) + 
     -           ((Alfa*MW*(-1 + (4 - 8*SW2)*SW2))/(4.D0*CW2**2*Pi*SW2)*
     -           (EL*Bval(bb00,iintL4(Gen4)))/SW + 
     -           (Alfa*MW*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -           (q12*EL*Bval(bb1,iintL4(Gen4)))/SW + 
     -           (Alfa*MW*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -           (EL*cintL1(Gen4))/SW + (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -           (EL*Bval(bb0,iintL4(Gen4))*Mf2(2,Gen4))/SW)/(-MZ2 +q12)
         end do
      end if

** fermion formfactors
      a1selfZfm = loopa1 + ctFM


** Last, and most annoyingly, the bosons:
      if (sector .eq. 3) then
         loopa1 = loopa1 + ((Alfa*MW)/(6.D0*Pi*SW2)*(q22*EL)/SW + 
     -        Bval(bb0,iintB2)*((5*Alfa*MW)/(4.D0*Pi*SW2)*(q22*EL)/SW + 
     -        (Alfa*MW*MW2*(CW2**2 - SW2**2))/(2.D0*CW2**2*Pi*SW2)*
     -        EL/SW) + (Alfa*MW*(-SW2**2 + CW2*(-9*CW2 + 2*SW2)))/
     -        (8.D0*CW2**2*Pi*SW2)*(cintB1*EL)/SW - 
     -        (Alfa*MW)/(16.D0*CW2**2*Pi*SW2)*((cintB2 +cintB3)*EL)/SW - 
     -        (Alfa*MW*MW2)/(4.D0*CW2**3*Pi*SW2)*
     -        (EL*Bval(bb0,iintB1))/SW + (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb00,iintB1))/SW + 
     -        (Alfa*MW*(CW2*(9*CW2 - 2*SW2) + SW2**2))/
     -        (4.D0*CW2**2*Pi*SW2)*(EL*Bval(bb00,iintB2))/SW + 
     -        (Alfa*MW)/(2.D0*Pi*SW2)*(q22*EL*Bval(bb1,iintB2))/SW)/
     -        (-MZ2 + q22) + ((Alfa*MW)/(6.D0*Pi*SW2)*(q12*EL)/SW + 
     -        Bval(bb0,iintB4)*((5*Alfa*MW)/(4.D0*Pi*SW2)*(q12*EL)/SW + 
     -        (Alfa*MW*MW2*(CW2**2 - SW2**2))/(2.D0*CW2**2*Pi*SW2)*
     -        EL/SW) + (Alfa*MW*(-SW2**2 + CW2*(-9*CW2 + 2*SW2)))/
     -        (8.D0*CW2**2*Pi*SW2)*(cintB1*EL)/SW - 
     -        (Alfa*MW)/(16.D0*CW2**2*Pi*SW2)*((cintB2 + cintB3)*EL)/SW - 
     -        (Alfa*MW*MW2)/(4.D0*CW2**3*Pi*SW2)*
     -        (EL*Bval(bb0,iintB3))/SW + (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -        (EL*Bval(bb00,iintB3))/SW + 
     -        (Alfa*MW*(CW2*(9*CW2 - 2*SW2) + SW2**2))/
     -        (4.D0*CW2**2*Pi*SW2)*(EL*Bval(bb00,iintB4))/SW + (Alfa*MW)/
     -        (2.D0*Pi*SW2)*(q12*EL*Bval(bb1,iintB4))/SW)/(-MZ2 + q12)
      end if
      
      
      if (debug .ne. 0) then
         write(*,*)'ZZ self energy loopa1=', loopa1
      end if


** Now we want the form factors themselves
      a1selfZ = a1selfZ + loopa1 + ct 

      if (debug .ne. 0) then
         write(*,*)'ZZ self energy a1 =', a1selfZ
         write(*,*)'ZZ FERMION contrn to self energy a1 =', a1selfZ
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the AZ self energy

      subroutine az_selfSM(a1selfAZ,q12,q22,q32,q1q2,debug,fm)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** These are the contributions from the counter-term (g-mu-nu) and the loop
** (g-mu-nu).
      double complex ct, loopa1

** Integral arguments
      double complex cint1(3), cint2(3)
      integer*8 iint1(3), iint2(3)
      double complex cintL1(3)
      integer*8 iintL1(3)
      double complex cintB1
      integer*8 iintB1


** This is the contribution to the formfactor a1 (g-mu-nu)
      double complex a1selfAZ


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions.  'fm' = 1 for AZ and 2 for AZ2
      integer fm
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


** Particle do loop parameters
      integer g4ini, Gen4

      
      double precision Alfa, Alfa2

      integer debug


** This is needed for looptools. 
#include "looptools.h"

** Parameter declarations
#include "mssm.inc"
#include "global.inc"


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions and 3 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Clearing the loop and counterterm values
      loopa1 = (0D0,0D0)  
      ct = (0D0,0D0)    


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
        ct = 1/(2.D0*CW2)*(q12*(-(MW*(dZAZ1*EL)) - MW*(dZZA1*EL)) + 
     -       MW*MZ2*(dZZA1*EL))/((-MZ2 + q12)*SW)

** counterterm for fermions
        ctFM = 1/(2.D0*CW2)*(q12*(-(MW*(dZAZ1fm*EL)) - MW*(dZZA1fm*EL)) + 
     -       MW*MZ2*(dZZA1fm*EL))/((-MZ2 + q12)*SW)
      
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'AZ self energy ct=', ct
         write(*,*)'AZ self energy fermion ct=', ctFM
      end if


** Now the loop

** Now for the arguments of the integrals.
      do Gen4 = g4ini, 3
         cint1(Gen4) = A0(Mf2(3,Gen4))
         cint2(Gen4) = A0(Mf2(4,Gen4))
         iint1(Gen4) = Bget(q12,Mf2(3,Gen4),Mf2(3,Gen4))
         iint2(Gen4) = Bget(q12,Mf2(4,Gen4),Mf2(4,Gen4))
      enddo

      if (sector .ne. 1) then
         do Gen4 = 1, 3
            cintL1(Gen4) = A0(Mf2(2,Gen4))
            iintL1(Gen4) = Bget(q12,Mf2(2,Gen4),Mf2(2,Gen4))
         enddo
      end if

      if (sector .eq. 3) then   
         cintB1 = A0(MW2)
         iintB1 = Bget(q12,MW2,MW2)
      end if

      
** First, the quarks:
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + 1/(12.D0*CW*CW2*Pi*SW2)*(Alfa*MW*(q12*EL*(-6*
     -        Bval(bb1,iint1(Gen4)) - 3*Bval(bb1,iint2(Gen4)) + SW2*(16*
     -        Bval(bb1,iint1(Gen4)) + 4*Bval(bb1,iint2(Gen4))))) + Alfa*
     -        MW*(EL*(12*Bval(bb00,iint1(Gen4)) + 6*
     -        Bval(bb00,iint2(Gen4)) - 6*cint1(Gen4) - 3*cint2(Gen4) + 
     -        SW2*(-32*Bval(bb00,iint1(Gen4))-8*Bval(bb00,iint2(Gen4)) + 
     -        16*cint1(Gen4) + 4*cint2(Gen4)))))/(-MZ2 + q12)
      end do


** Now, the leptons ...
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            loopa1 = loopa1 + 
     -           1/(4.D0*CW*CW2*Pi*SW2)*(Alfa*MW*(q12*EL*
     -           (-Bval(bb1,iintL1(Gen4)) + 4*SW2*
     -           Bval(bb1,iintL1(Gen4)))) + Alfa*MW*(EL*
     -           (2*Bval(bb00,iintL1(Gen4)) - cintL1(Gen4) + SW2*(-8*
     -           Bval(bb00,iintL1(Gen4))+ 4*cintL1(Gen4)))))/(-MZ2 +q12)
         end do
      end if


** fermion formfactor
      if (fm .eq. 1) then
         a1selfAZfm = loopa1 + ctFM 
      else
         a1selfAZ2fm = loopa1 + ctFM 
      end if


** Last, and most annoyingly, the bosons:
      if (sector .eq. 3) then
         loopa1 = loopa1+(Bval(bb0,iintB1)*(-((5*Alfa*CW*MW)/(4.D0*CW2*
     -        Pi*SW2)*(q12*EL))-(Alfa*MW*MW2)/(2.D0*CW*CW2*Pi*SW2)*EL) - 
     -        (Alfa*CW*MW)/(2.D0*CW2*Pi*SW2)*(q12*EL*Bval(bb1,iintB1)) + 
     -        Alfa/(4.D0*CW*CW2*Pi*SW2)*(cintB1*(5*CW2*MW*EL - MW*SW2*
     -        EL)) + Alfa/(2.D0*CW*CW2*Pi*SW2)*(Bval(bb00,iintB1)*(-(5*
     -        CW2*MW*EL) + MW*SW2*EL)))/(-MZ2 + q12)
      end if
      
      
      if (debug .ne. 0) then
         write(*,*)'AZ self energy loopa1=', loopa1
      end if


** Now we want the form factors themselves
      a1selfAZ = loopa1 + ct 
      
      
      if (debug .ne. 0) then
         write(*,*)'AZ self energy a1 =', a1selfAZ
         write(*,*)'AZ FERMION self energy a1 =', a1selfAZfm, a1selfAZ2fm
      end if


      end


*******************************************************************************


*************************************************
***                   SUSY                    ***
*************************************************


******************************************************************************* 

*** This subroutine calculates the tree-level part of the a1 form factor for 
*** the T-mu-nu anomalous coupling for the VV-higgs vertex in the MSSM

      subroutine treeSUSY(treeW,treeZ,debug,higgno)

      implicit none


** Value of tree level a1 coupling for W and Z bosons coupling.  The first 
** entry in each array is the light higgs, the second entry is the heavy higgs.
      double complex treeW(3), treeZ(3)


** This parameter defines which of the neutral higgs bosons is being considered.
** It is set in vbfnlo.dat using HMASS
      integer higgno

      integer debug


** Parameter declarations
#include "mssm.inc"


** The light higgs
      treeW(1) = MW*(EL*SBA)/SW
      treeZ(1) = 1/CW2*(MW*(EL*SBA)/SW) 

** The heavy higgs    
      treeW(2) = MW*(CBA*EL)/SW
      treeZ(2) = 1/CW2*(MW*(EL*CBA)/SW) 

** The cp odd higgs coupling is 0 at tree level
      treeW(3) = (0D0,0D0)
      treeZ(3) = (0D0,0D0)


      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'WWh tree=', treeW(1)
         write(*,*)'ZZh tree=', treeZ(1)
         write(*,*)'WWH tree=', treeW(2)
         write(*,*)'ZZH tree=', treeZ(2)
      end if

      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the WWh/H/A vertex 
*** in susy.

      subroutine ww_vtxSUSY(a1vtxW,a2vtxW,a3vtxW,q12,q22,q32,q1q2,debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the z bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Integral arguments (suffix Q -> 'quark', L -> 'leptonic', B -> 'bosonic')
      integer*8 iintQ1(3), iintQ2(3), iintQ3(3), iintQ4(3), iintQ5(3,2,2)
      integer*8 iintQ6(3,2,2), iintQ7(3,2,2,2), iintQ8(3,2,2,2)
      integer*8 iintL1(3), iintL2(3), iintL3(3), iintL4(3,2) 
      integer*8 iintL5(3,2,2), iintL6(3,2,2)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6, iintB7 
      integer*8 iintB8, iintB9, iintB10, iintB11(3), iintB12(3) 
      integer*8 iintB13(3), iintB14(3), iintB15(3), iintB16(3), iintB17(3)
      integer*8 iintB18(3), iintB19(3,3), iintB20(3,3), iintB21(3,3)
      integer*8 iintB22(2,4,4), iintB23(4,4), iintB24(2,2,4), iintB25(2,2)

** These are the contributions from the counter-term (g-mu-nu) and the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, cth0, ctHH, loopa1, loopa2, loopa3


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for WW-higgs.  The first entry in each array is the light higgs, the second 
** entry is the heavy higgs.
      double complex a1vtxW, a2vtxW, a3vtxW


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions
      double complex ctFM, cth0FM, ctHHfm
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


** Particle do loop parameters
      integer g4ini, Gen4, Sfe5, Sfe4, Sfe6
      integer Hig4, Hig5, Cha4, Cha5, Neu4, Neu5

** Dummy parameters (q for 'quark', l for 'lepton', b for 'boson')
      double complex qup1, qup2, qup3, qup4, qup5, qup6, qup7
      double complex qup8, qup9, qup10, qup11, qup12, qup13, qup14
      double complex qup15, qup16, qup17, qup18, qup19, qup20, qup21
      double complex qup22, qup23, qup24, qup25, qup26, qup27, qup28
      double complex lup1, lup2, lup3, lup4, lup5, lup6, lup7
      double complex lup8, lup9, lup10, lup11, lup12, lup13, lup14
      double complex bup1, bup2, bup3, bup4, bup5, bup6, bup7
      double complex bup8, bup9, bup10, bup11, bup12, bup13, bup14
      double complex bup15, bup16, bup17, bup18, bup19, bup20, bup21
      double complex bup22, bup23, bup24, bup25, bup26, bup27, bup28
      double complex bup29, bup30, bup31, bup32, bup33, bup34, bup35
      double complex bup36, bup37, bup38, bup39, bup40, bup41, bup42
      double complex bup43, bup44, bup45, bup46, bup47, bup48, bup49
      double complex bup50, bup51, bup52


      double precision Alfa, Alfa2

** Higgs mixing matrices
      double complex UHggs(3,3), ZHggs(3,3)

      integer debug

** This is needed for looptools. 
#include "looptools.h"

** Parameter declarations
#include "mssm.inc"
#include "global.inc"


** UHggs and ZHggs: set to 1 for debugging
      do Hig4 = 1, 3
         do Hig5 = 1, 3
            if (Hig4 .eq. Hig5) then
               UHggs(Hig4,Hig5) = (1D0,0D0)
*               ZHggs(Hig4,Hig5) = (1D0,0D0)
            else
               UHggs(Hig4,Hig5) = (0D0,0D0)
*               ZHggs(Hig4,Hig5) = (0D0,0D0)
            end if
         end do
      end do


      if (propLEVEL .ne. 2) then
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               if (Hig4 .eq. Hig5) then
                  ZHggs(Hig4,Hig5) = (1D0,0D0)
               else
                  ZHggs(Hig4,Hig5) = (0D0,0D0)
               end if
            end do
         end do
      else
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               ZHggs(Hig4,Hig5) = Hgs(Hig4,Hig5)
            end do
         end do
      end if


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Clearing the values of the formfactors
      a1vtxW = (0D0,0D0)
      a2vtxW = (0D0,0D0)
      a3vtxW = (0D0,0D0)


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions, 3 for all SM-type and
** sfermions, and 4, 5, 6 for all particles included in the hVV vertex.
      if (sector .eq. 1) then
         g4ini = 3
      else 
         g4ini = 1
      end if


** Now for the arguments of the integrals. 
** First the quarks and squarks
      do Gen4 = g4ini, 3
         iintQ1(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ2(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),
     -        Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ3(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
         iintQ4(Gen4) = Cget(q32,q12,q22,Mf2(3,Gen4),
     -        Mf2(3,Gen4),Mf2(4,Gen4))
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3      
               iintQ5(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,4,Gen4),MSf2(Sfe5,4,Gen4))
               iintQ6(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,3,Gen4)) 
            enddo
         enddo
      enddo

      do Sfe6 = 1,2
         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = g4ini, 3
                  iintQ7(Gen4,Sfe4,Sfe5,Sfe6) = Cget(q12,q32,q22,
     -                 MSf2(Sfe4,3,Gen4),MSf2(Sfe5,4,Gen4),
     -                 MSf2(Sfe6,4,Gen4))
                  iintQ8(Gen4,Sfe4,Sfe5,Sfe6) = Cget(q12,q32,q22,
     -                 MSf2(Sfe4,4,Gen4),MSf2(Sfe5,3,Gen4),
     -                 MSf2(Sfe6,3,Gen4))
               enddo
            enddo
         enddo
      enddo

** Now the leptonic and sleptonic integrals
      if (sector .ne. 1) then
         do Gen4 = 1,3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL2(Gen4) = Cget(q22,q32,q12,0.D0,Mf2(2,Gen4),
     -           Mf2(2,Gen4))
            iintL3(Gen4) = Bget(q32,MSf2(1,1,Gen4),MSf2(1,1,Gen4))
         enddo

         do Sfe4 = 1,2
            do Gen4 = 1,3 
               iintL4(Gen4,Sfe4) = 
     -              Cget(q32,q12,q22,MSf2(1,1,Gen4),MSf2(1,1,Gen4),
     -              MSf2(Sfe4,2,Gen4))              
            enddo
         enddo

         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  iintL5(Gen4,Sfe4,Sfe5) = 
     -                 Bget(q32,MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4))
                  iintL6(Gen4,Sfe4,Sfe5) = Cget(q12,q32,q22,
     -                 MSf2(1,1,Gen4),MSf2(Sfe4,2,Gen4),
     -                 MSf2(Sfe5,2,Gen4))
               enddo
            enddo
         enddo
      end if

** The 'bosonic' integrals (i.e. loops involving sm bosons and their 
** superpartners)...
      if (sector .ge. 3) then

** SM-type integrals
         iintB1 = Bget(q12,0.D0,MW2)
         iintB2 = Bget(q22,0.D0,MW2)
         iintB3 = Bget(q12,MW2,MZ2)
         iintB4 = Bget(q22,MW2,MZ2)
         iintB5 = Cget(q22,q32,q12,0.D0,MW2,MW2)
         iintB6 = Bget(q32,MHp2,MHp2)
         iintB7 = Bget(q32,MW2,MW2)
         iintB8 = Cget(q22,q32,q12,MW2,MZ2,MZ2)
         iintB9 = Cget(q32,q12,q22,MW2,MW2,MZ2)
         iintB10 = Bget(q32,MZ2,MZ2)

** neutral higgs integrals ..
         do Hig4 = 1,3
            iintB11(Hig4) = Cget(q32,q22,q12,MHp2,MW2,MHiggs2(Hig4))
            iintB12(Hig4) = Cget(q32,q12,q22,MHp2,MW2,MHiggs2(Hig4))
            iintB13(Hig4) = Cget(q32,q12,q22,MW2,MW2,MHiggs2(Hig4))
            iintB14(Hig4) = Cget(q32,q12,q22,MHp2,MHp2,MHiggs2(Hig4))
            iintB15(Hig4) = Cget(q12,q32,q22,MW2,MZ2,MHiggs2(Hig4))
            iintB16(Hig4) = Cget(q22,q32,q12,MW2,MZ2,MHiggs2(Hig4))
            iintB17(Hig4) = Bget(q12,MW2,MHiggs2(Hig4))
            iintB18(Hig4) = Bget(q22,MW2,MHiggs2(Hig4))
         enddo

         do Hig5 = 1,3
            do Hig4 = 1,3
               iintB19(Hig4,Hig5) = 
     -              Cget(q12,q32,q22,MW2,MHiggs2(Hig4),MHiggs2(Hig5))
               iintB20(Hig4,Hig5) = 
     -              Cget(q12,q32,q22,MHp2,MHiggs2(Hig4),MHiggs2(Hig5))
               iintB21(Hig4,Hig5) = Bget(q32,MHiggs2(Hig4),MHiggs2(Hig5))
            enddo
         enddo

** chargino and neutralino integrals
         if (sector .ge. 4) then
            do Cha5 = 1,2
               do Cha4 = 1,2
                  iintB25(Cha4,Cha5) = Bget(q32,MCha2(Cha4),MCha2(Cha5))
               enddo
            enddo

            do Neu5 = 1,4
               do Neu4 = 1,4
                  iintB23(Neu4,Neu5) = Bget(q32,MNeu2(Neu4),MNeu2(Neu5))
               enddo
            enddo

            do Neu4 = 1,4
               do Cha5 = 1,2
                  do Cha4 = 1,2
                     iintB24(Cha4,Cha5,Neu4) = Cget(q32,q22,q12,
     -                    MCha2(Cha4),MCha2(Cha5),MNeu2(Neu4))
                  enddo
               enddo
            enddo

            do Neu5 = 1,4
               do Neu4 = 1,4
                  do Cha4 = 1,2
                     iintB22(Cha4,Neu4,Neu5) = Cget(q12,q32,q22,
     -                    MCha2(Cha4),MNeu2(Neu4),MNeu2(Neu5))
                  enddo
               enddo
            enddo
         end if

      end if


** Clearing the loop values
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0) 
      ct = (0D0,0D0) 
      cth0 = (0D0,0D0) 
      ctHH = (0D0,0D0) 


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterms
      cth0 = MW/(2.D0*MW2*SW2)*(EL*(MW2*SW2*(CBA*dZh0HH) + 2*MW2*SW2*
     -     (CB*CBA*dTB1*SB) + SW2*(dMWsq1*SBA) + 2*MW2*SW2*(dZet1*SBA) + 
     -     MW2*SW2*(dZh0*SBA) + 2*MW2*SW2*(dZW1*SBA) + 
     -     2*MW2*(dSW1*SBA*SW)))/SW

      ctHH = MW/(2.D0*MW2*SW2)*(EL*(MW2*SW2*(dZh0HH*SBA) - 2*MW2*SW2*
     -     (CB*dTB1*SB*SBA) + CBA*(SW2*dMWsq1 + 2*MW2*SW2*dZet1 + 
     -     MW2*SW2*dZHH + 2*MW2*SW2*dZW1 + 2*MW2*(dSW1*SW))))/SW

      if (sector .eq. 4) then
         cth0 = MW/(2.D0*MW2*SW2)*(EL*(MW2*SW2*(CBA*dZh0HHfull) + 2*MW2*
     -        SW2*(CB*CBA*dTB1full*SB) + SW2*(dMWsq1full*SBA) + 2*MW2*
     -        SW2*(dZet1full*SBA) + MW2*SW2*(dZh0full*SBA) + 2*MW2*SW2*
     -        (dZW1full*SBA) + 2*MW2*(dSW1full*SBA*SW)))/SW
         
         ctHH = MW/(2.D0*MW2*SW2)*(EL*(MW2*SW2*(dZh0HHfull*SBA) - 2*MW2*
     -        SW2*(CB*dTB1full*SB*SBA) + CBA*(SW2*dMWsq1full + 2*MW2*
     -        SW2*dZet1full + MW2*SW2*dZHHfull + 2*MW2*SW2*dZW1full + 2*
     -        MW2*(dSW1full*SW))))/SW
      end if
      
      ct = cth0*ZHggs(HiggsType,1) + ctHH*ZHggs(HiggsType,2)

** for the fermions
      cth0FM = MW/(2.D0*MW2*SW2)*(EL*(MW2*SW2*(CBA*dZh0HHfm) + 2*MW2*SW2*
     -     (CB*CBA*dTB1fm*SB) + SW2*(dMWsq1fm*SBA) + 2*MW2*SW2*(dZet1fm*SBA) + 
     -     MW2*SW2*(dZh0fm*SBA) + 2*MW2*SW2*(dZW1fm*SBA) + 
     -     2*MW2*(dSW1fm*SBA*SW)))/SW

      ctHHfm = MW/(2.D0*MW2*SW2)*(EL*(MW2*SW2*(dZh0HHfm*SBA) - 2*MW2*SW2*
     -     (CB*dTB1fm*SB*SBA) + CBA*(SW2*dMWsq1fm + 2*MW2*SW2*dZet1fm + 
     -     MW2*SW2*dZHHfm + 2*MW2*SW2*dZW1fm + 2*MW2*(dSW1fm*SW))))/SW
      
      ctFM = cth0FM*ZHggs(HiggsType,1) + ctHHfm*ZHggs(HiggsType,2)

      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'WW-higgs ct=', ct
         write(*,*)'WWh0 =', cth0
         write(*,*)'WWHH =', ctHH
         write(*,*)'WW-higgs FERMION contn to ct=', ctFM
         write(*,*)'FERMION contn to WWh0 ct =', cth0FM
         write(*,*)'FERMION contn to WWHH ct=', ctHHfm
      end if


** Now the loop - first the quarks
      do Gen4 = g4ini, 3
         qup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
         qup2 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
         qup3 = -(cI*(SB*ZHggs(HiggsType,3)))+ SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)
         qup4 = -(cI*(SB*ZHggs(HiggsType,3))) + CA*TB*
     -        ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)
         qup5 = cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)
         qup6 = cI*(SB*ZHggs(HiggsType,3)) +CA*TB*ZHggs(HiggsType,1) + 
     -        SA*TB*ZHggs(HiggsType,2)

         loopa1 = loopa1 + 1/(16.D0*MW*Pi*SW2)*(Cval(cc0,iintQ4(Gen4))*
     -        (-(3*Alfa*((q12+q22-q32)*qup4*EL*Mf2(3,Gen4))) - 6*Alfa*
     -        (q22*qup6*EL*Mf2(3,Gen4)) + 12*Alfa*(qup2*EL*TB*
     -        Mf2(3,Gen4)*Mf2(4,Gen4))))/(SB*SW*TB) + (3*Alfa)/(4.D0*MW*
     -        Pi*SW2)*(qup2*EL*Bval(bb0,iintQ3(Gen4))*Mf2(3,Gen4))/(SB*
     -        SW) - (3*Alfa)/(2.D0*MW*Pi*SW2)*(qup2*EL*
     -        Cval(cc00,iintQ4(Gen4))*Mf2(3,Gen4))/(SB*SW) + (3*Alfa)/
     -        (16.D0*MW*Pi*SW2)*(((q12-q22+q32)*qup4 + (q12-q22-q32)*
     -        qup6)*EL*Cval(cc1,iintQ4(Gen4))*Mf2(3,Gen4))/(SB*SW*TB) - 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*(((q12+q22-q32)*qup4 + 2*q22*
     -        qup6)*EL*Cval(cc2,iintQ4(Gen4))*Mf2(3,Gen4))/(SB*SW*TB) - 
     -        (3*Alfa)/(4.D0*MW*Pi*SW2)*(qup1*EL*Bval(bb0,iintQ1(Gen4))*
     -        Mf2(4,Gen4))/(CB*SW) + (3*Alfa)/(2.D0*MW*Pi*SW2)*(qup1*EL*
     -        Cval(cc00,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) - (3*Alfa)/
     -        (16.D0*MW*Pi*SW2)*(((q12+q22-q32)*qup3 + 2*q22*qup5)*EL*
     -        Cval(cc1,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) - (3*Alfa)/
     -        (16.D0*MW*Pi*SW2)*((2*q12*qup3 + (q12+q22-q32)*qup5)*EL*
     -        Cval(cc2,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) - (3*Alfa)/
     -        (4.D0*MW*Pi*SW2)*(qup1*EL*Cval(cc0,iintQ2(Gen4))*
     -        Mf2(3,Gen4)*Mf2(4,Gen4))/(CB*SW)


         loopa2 = loopa2 + 1/(16.D0*MW*Pi*SW2)*
     -        (Cval(cc0,iintQ4(Gen4))*(-(3*Alfa*(2*qup4*EL*
     -        Mf2(3,Gen4)))))/(SB*SW*TB) + (3*Alfa)/(16.D0*MW*Pi*SW2)*
     -        (((-6)*qup4 + (-2)*qup6)*EL*Cval(cc1,iintQ4(Gen4))*
     -        Mf2(3,Gen4))/(SB*SW*TB) - (3*Alfa)/(2.D0*MW*Pi*SW2)*
     -        (qup2*EL*(Cval(cc11,iintQ4(Gen4)) + 
     -        Cval(cc12,iintQ4(Gen4)))*Mf2(3,Gen4))/(SB*SW) - (3*Alfa)/
     -        (16.D0*MW*Pi*SW2)*((2*qup4)*EL*Cval(cc2,iintQ4(Gen4))*
     -        Mf2(3,Gen4))/(SB*SW*TB) - (3*Alfa)/(16.D0*MW*Pi*SW2)*
     -        ((2*qup3)*EL*Cval(cc1,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) - 
     -        (3*Alfa)/(2.D0*MW*Pi*SW2)*(qup1*EL*
     -        Cval(cc12,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) - (3*Alfa)/
     -        (16.D0*MW*Pi*SW2)*((2*qup5)*EL*
     -        Cval(cc2,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) 


         loopa3 = loopa3 + 1/(16.D0*MW*Pi*SW2)*(Cval(cc0,iintQ4(Gen4))*
     -        (-(3*Alfa*((-2)*qup4*EL*Mf2(3,Gen4)))))/(SB*SW*TB) + (3*
     -        Alfa)/(16.D0*MW*Pi*SW2)*((2*qup4 + 2*qup6)*EL*
     -        Cval(cc1,iintQ4(Gen4))*Mf2(3,Gen4))/(SB*SW*TB) - (3*Alfa)/
     -        (16.D0*MW*Pi*SW2)*(((-2)*qup4)*EL*Cval(cc2,iintQ4(Gen4))*
     -        Mf2(3,Gen4))/(SB*SW*TB)  - (3*Alfa)/(16.D0*MW*Pi*SW2)*
     -        ((2*qup3)*EL*Cval(cc1,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) - 
     -        (3*Alfa)/(16.D0*MW*Pi*SW2)*(((-2)*qup5)*EL*
     -        Cval(cc2,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) 
      end do

         
** Contributions from squarks         
      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3
               qup7 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
               qup8 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
               qup9 = -(SAB*ZHggs(HiggsType,1))+CAB*ZHggs(HiggsType,2)
               qup10 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
               qup11 = -(cI*(SB*ZHggs(HiggsType,3))) +
     -              SA*ZHggs(HiggsType,1)- CA*ZHggs(HiggsType,2)
               qup12 = cI*(SB*ZHggs(HiggsType,3)) + 
     -              SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)

               loopa1 = loopa1 - Alfa/(16.D0*CW*MW*Pi*SW2)*(EL*
     -              Bval(bb0,iintQ6(Gen4,Sfe4,Sfe5))*USf(Sfe4,1,3,Gen4)*
     -              (USf(Sfe5,1,3,Gen4)*(3*CW*(Mf(3,Gen4)*
     -              USfC(Sfe4,2,3,Gen4)*(qup11*MUE*TB +AfC(3,Gen4,Gen4)*
     -              (-(cI*(SB*ZHggs(HiggsType,3))) + CA*TB*
     -              ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)))) + 
     -              TB*(MW*MZ*(-3 + 4*SW2)*
     -              (qup8*SB) + 6*CW*(qup10*Mf2(3,Gen4)))*
     -              USfC(Sfe4,1,3,Gen4)) + USf(Sfe5,2,3,Gen4)*(3*CW*
     -              (Mf(3,Gen4)*USfC(Sfe4,1,3,Gen4)*(qup12*MUEC*TB + 
     -              Af(3,Gen4,Gen4)*(cI*(SB*ZHggs(HiggsType,3)) +CA*TB*
     -              ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)))) + 
     -              2*TB*(2*MW*MZ*SW2*(qup9*SB) + 3*CW*(qup10*
     -              Mf2(3,Gen4)))*USfC(Sfe4,2,3,Gen4)))*
     -              USfC(Sfe5,1,3,Gen4))/(SB*SW*TB) + Alfa/(16.D0*CW*MW*
     -              Pi*SW2)*(EL*Bval(bb0,iintQ5(Gen4,Sfe4,Sfe5))*
     -              USf(Sfe4,1,4,Gen4)*(USf(Sfe5,1,4,Gen4)*(3*CW*
     -              (Mf(4,Gen4)*USfC(Sfe4,2,4,Gen4)*(qup12*
     -              AfC(4,Gen4,Gen4) +MUE*(cI*(CB*ZHggs(HiggsType,3)) + 
     -              CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)))) + 
     -              (MW*MZ*(-3 + 2*SW2)*(CB*qup8) + 6*CW*(qup7*
     -              Mf2(4,Gen4)))*USfC(Sfe4,1,4,Gen4))+
     -              USf(Sfe5,2,4,Gen4)*(3*CW*(Mf(4,Gen4)*
     -              USfC(Sfe4,1,4,Gen4)*(qup11*Af(4,Gen4,Gen4) + MUEC*
     -              (-(cI*(CB*ZHggs(HiggsType,3))) + CA*
     -              ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)))) + 2*
     -              (MW*MZ*SW2*(CB*qup9) + 3*CW*(qup7*Mf2(4,Gen4)))*
     -              USfC(Sfe4,2,4,Gen4)))*USfC(Sfe5,1,4,Gen4))/(CB*SW)
            end do
         end do
      end do


** More contributions from squarks
      do Sfe6 = 1,2
         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = g4ini, 3 
                  qup13 = SA*ZHggs(HiggsType,1) - CA*
     -                 ZHggs(HiggsType,2)
                  qup14 = SAB*ZHggs(HiggsType,1) - CAB*
     -                 ZHggs(HiggsType,2)
                  qup15 = -(SAB*ZHggs(HiggsType,1)) + CAB*
     -                 ZHggs(HiggsType,2)
                  qup16 = CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                  qup17 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -                 CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                  qup18 = cI*(CB*ZHggs(HiggsType,3)) +
     -                 CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                  qup19 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -                 SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                  qup20 = -(cI*(SB*ZHggs(HiggsType,3))) + CA*TB*
     -                 ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)
                  qup21 = cI*(SB*ZHggs(HiggsType,3)) +
     -                 SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                  qup22 = cI*(SB*ZHggs(HiggsType,3)) + CA*TB*
     -                 ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)
                  qup23 = MW*MZ*(-3 + 4*SW2)*(qup14*SB) + 
     -                 6*CW*(qup16*Mf2(3,Gen4))
                  qup24 = MW*MZ*(-3 + 2*SW2)*(CB*qup14) + 
     -                 6*CW*(qup13*Mf2(4,Gen4))
                  qup25 = 3*CW*((qup18*MUE + qup21*AfC(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe5,2,4,Gen4)) + 
     -                 qup24*USfC(Sfe5,1,4,Gen4)
                  qup26 = 3*CW*((qup17*MUEC + qup19*Af(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe5,1,4,Gen4)) + 2*(MW*MZ*SW2*
     -                 (CB*qup15) + 3*CW*(qup13*Mf2(4,Gen4)))*
     -                 USfC(Sfe5,2,4,Gen4)
                  qup27 = 3*CW*((qup19*MUE*TB + qup20*AfC(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe6,2,3,Gen4)) + 
     -                 qup23*TB*USfC(Sfe6,1,3,Gen4)
                  qup28 = 3*CW*((qup21*MUEC*TB + qup22*Af(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe6,1,3,Gen4)) + 
     -                 2*TB*(2*MW*MZ*SW2*(qup15*SB) + 3*CW*
     -                 (qup16*Mf2(3,Gen4)))*USfC(Sfe6,2,3,Gen4)

                  loopa1 = loopa1 + Alfa/(4.D0*CW*MW*Pi*SW2)*
     -                 (EL*Cval(cc00,iintQ8(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 USf(Sfe4,1,4,Gen4)*(qup27*USf(Sfe5,1,3,Gen4) + 
     -                 qup28*USf(Sfe5,2,3,Gen4))*
     -                 USf(Sfe6,1,3,Gen4)*USfC(Sfe4,1,4,Gen4)*
     -                 USfC(Sfe5,1,3,Gen4))/(SB*SW*TB) - 
     -                 Alfa/(4.D0*CW*MW*Pi*SW2)*
     -                 (EL*Cval(cc00,iintQ7(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 USf(Sfe4,1,3,Gen4)*USf(Sfe5,1,4,Gen4)*
     -                 (qup25*USf(Sfe6,1,4,Gen4) + qup26*
     -                 USf(Sfe6,2,4,Gen4))*USfC(Sfe4,1,3,Gen4)*
     -                 USfC(Sfe6,1,4,Gen4))/(CB*SW) 
 

                  loopa2 = loopa2 - Alfa/(4.D0*CW*MW*Pi*SW2)*
     -                 (EL*Cval(cc12,iintQ8(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 USf(Sfe4,1,4,Gen4)*(qup27*USf(Sfe5,1,3,Gen4) + 
     -                 qup28*USf(Sfe5,2,3,Gen4))*
     -                 USf(Sfe6,1,3,Gen4)*USfC(Sfe4,1,4,Gen4)*
     -                 USfC(Sfe5,1,3,Gen4))/(SB*SW*TB) + 
     -                 Alfa/(4.D0*CW*MW*Pi*SW2)*
     -                 (EL*Cval(cc12,iintQ7(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 USf(Sfe4,1,3,Gen4)*USf(Sfe5,1,4,Gen4)*
     -                 (qup25*USf(Sfe6,1,4,Gen4) + qup26*
     -                 USf(Sfe6,2,4,Gen4))*
     -                 USfC(Sfe4,1,3,Gen4)*USfC(Sfe6,1,4,Gen4))/(CB*SW)
               enddo
            enddo
         enddo
      enddo


** Now the leptons and sleptons
      if (sector .ne. 1) then

** Firstly: leptons
         do Gen4 = 1,3
            lup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
            lup2 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
            lup3 = cI*(SB*ZHggs(HiggsType,3)) +SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)

            loopa1 = loopa1 - Alfa/(4.D0*MW*Pi*SW2)*(lup1*EL*
     -           Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4))/(CB*SW) + Alfa/
     -           (2.D0*MW*Pi*SW2)*(lup1*EL*Cval(cc00,iintL2(Gen4))*
     -           Mf2(2,Gen4))/(CB*SW) - Alfa/(16.D0*MW*Pi*SW2)*
     -           (((q22+q12-q32)*lup2 + 2*q22*lup3)*EL*
     -           Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) - Alfa/
     -           (16.D0*MW*Pi*SW2)*((2*q12*lup2 + (q22+q12-q32)*lup3)*
     -           EL*Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) + 
     -           (Alfa*MZ)/(16.D0*CW*Pi*SW2)*(EL*Bval(bb0,iintL3(Gen4))*
     -           (SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)))/SW

            loopa2 = loopa2 - Alfa/(16.D0*MW*Pi*SW2)*((2*lup2)*EL*
     -           Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) - Alfa/
     -           (2.D0*MW*Pi*SW2)*(lup1*EL*Cval(cc12,iintL2(Gen4))*
     -           Mf2(2,Gen4))/(CB*SW) - Alfa/(16.D0*MW*Pi*SW2)*
     -           ((2*lup3)*EL*Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4))/
     -           (CB*SW) 

            loopa3 = loopa3 - Alfa/(16.D0*MW*Pi*SW2)*((2*lup2)*EL*
     -           Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) - Alfa/
     -           (16.D0*MW*Pi*SW2)*(((-2)*lup3)*EL*
     -           Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) 
         enddo

         
** Now: sleptons
         do Sfe4 = 1,2
            do Gen4 = 1,3
               lup4 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)

               loopa1 = loopa1 - (Alfa*MZ)/(4.D0*CW*Pi*SW2)*
     -              (lup4*EL*Cval(cc00,iintL4(Gen4,Sfe4))*
     -              USf(Sfe4,1,2,Gen4)*USfC(Sfe4,1,2,Gen4))/SW 

               loopa2 = loopa2 - (Alfa*MZ)/(4.D0*CW*Pi*SW2)*
     -              (lup4*EL*(Cval(cc1,iintL4(Gen4,Sfe4)) + 
     -              Cval(cc11,iintL4(Gen4,Sfe4)) + 
     -              Cval(cc12,iintL4(Gen4,Sfe4)))*USf(Sfe4,1,2,Gen4)*
     -              USfC(Sfe4,1,2,Gen4))/SW
            enddo
         enddo
         
         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  lup5 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                  lup6 = SAB*ZHggs(HiggsType,1)-CAB*ZHggs(HiggsType,2)
                  lup7 = -(SAB*ZHggs(HiggsType,1)) + CAB*
     -                 ZHggs(HiggsType,2)
                  lup8 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -                 CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                  lup9 = cI*(CB*ZHggs(HiggsType,3)) + 
     -                 CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                  lup10 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -                 SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                  lup11 = cI*(SB*ZHggs(HiggsType,3)) + SA*
     -                 ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                  lup12 = MW*MZ*(-1 + 2*SW2)*(CB*lup6) + 2*CW*(lup5*
     -                 Mf2(2,Gen4))
                  lup13 = CW*((lup9*MUE + lup11*AfC(2,Gen4,Gen4))*
     -                 Mf(2,Gen4)*USfC(Sfe4,2,2,Gen4)) + 
     -                 lup12*USfC(Sfe4,1,2,Gen4)
                  lup14 = CW*((lup8*MUEC + lup10*Af(2,Gen4,Gen4))*
     -                 Mf(2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 
     -                 2*(MW*MZ*SW2*(CB*lup7) + CW*(lup5*Mf2(2,Gen4)))*
     -                 USfC(Sfe4,2,2,Gen4)

                  loopa1 = loopa1 + Alfa/(16.D0*CW*MW*Pi*SW2)*(EL*
     -                 Bval(bb0,iintL5(Gen4,Sfe4,Sfe5))*
     -                 USf(Sfe4,1,2,Gen4)*(lup13*USf(Sfe5,1,2,Gen4) + 
     -                 lup14*USf(Sfe5,2,2,Gen4))*USfC(Sfe5,1,2,Gen4))/
     -                 (CB*SW) - Alfa/(4.D0*CW*MW*Pi*SW2)*
     -                 (EL*Cval(cc00,iintL6(Gen4,Sfe4,Sfe5))*
     -                 USf(Sfe4,1,2,Gen4)*(lup13*USf(Sfe5,1,2,Gen4) + 
     -                 lup14*USf(Sfe5,2,2,Gen4))*USfC(Sfe5,1,2,Gen4))/
     -                 (CB*SW) 

                  loopa2 = loopa2 + Alfa/(4.D0*CW*MW*Pi*SW2)*
     -                 (EL*Cval(cc12,iintL6(Gen4,Sfe4,Sfe5))*
     -                 USf(Sfe4,1,2,Gen4)*(lup13*USf(Sfe5,1,2,Gen4) + 
     -                 lup14*USf(Sfe5,2,2,Gen4))*
     -                 USfC(Sfe5,1,2,Gen4))/(CB*SW)
               enddo
            enddo
         enddo

      end if


** fermion formfactors       
      a2vtxWfm = -loopa2
      a1vtxWfm = loopa1 + (loopa2*q1q2) + ctFM
      a3vtxWfm = loopa3*(0.D0,-1.D0)
      

** And finally the SM bosons and their superpartners
      if (sector .ge. 3) then
** SM-type
         bup1 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
         bup2 = SBA*ZHggs(HiggsType,1) + CBA*ZHggs(HiggsType,2)

         loopa1 = loopa1 + Bval(bb0,iintB7)*((Alfa*MW)/(16.D0*CW2*Pi*
     -        SW2)*(C2B*bup1*EL)/SW - (Alfa*MW*(3 - 2*CW2 - 2*SW2))/
     -        (4.D0*Pi*SW2)*(bup2*EL)/SW) - (Alfa*MW)/(8.D0*Pi)*
     -        (bup2*EL*(Bval(bb0,iintB1) + Bval(bb0,iintB2)))/SW - 
     -        (Alfa*MW)/(8.D0*CW2*Pi)*(bup2*EL*SW*(Bval(bb0,iintB3) + 
     -        Bval(bb0,iintB4))) + (2*Alfa*MW)/Pi*(bup2*EL*
     -        Cval(cc00,iintB5))/SW + 1/(8.D0*CW2*Pi*SW2)*
     -        (Cval(cc1,iintB8)*(Alfa*CW2*MW*((q12+3*q22-q32)*bup2*
     -        EL)))/SW + 1/(8.D0*CW2*Pi*SW2)*(Cval(cc2,iintB8)*(3*Alfa*
     -        CW2*MW*((q12+1/3.D0*q22-1/3.D0*q32)*bup2*EL)))/SW + 1/
     -        (8.D0*Pi)*(Cval(cc1,iintB5)*(Alfa*MW*((q12+3*q22-q32)*
     -        bup2*EL) + Alfa*MW*((q12-q22-q32)*bup2*EL) - Alfa*MW*
     -        ((q12-q22-q32)*bup2*EL)))/SW + 1/(8.D0*Pi*SW2)*
     -        (Cval(cc2,iintB9)*(-(Alfa*CW2*MW*((q12+3*q22-q32)*bup2*
     -        EL)) + Alfa*MW*SW2*((q12-q22-q32)*bup2*EL) - Alfa*MW*SW2*
     -        ((q12-q22-q32)*bup2*EL)))/SW + 1/(8.D0*Pi)*
     -        (Cval(cc2,iintB5)*(3*Alfa*MW*((q12+1/3.D0*q22-1/3.D0*q32)*
     -        bup2*EL) + Alfa*MW*((q12-q22+q32)*bup2*EL) - Alfa*MW*
     -        ((q12-q22+q32)*bup2*EL)))/SW + 1/(8.D0*CW2*Pi)*
     -        (Cval(cc0,iintB5)*(5*Alfa*CW2*MW*((q12+q22-q32)*bup2*EL) - 
     -        Alfa*CW2*MW*((q12-q32)*bup2*EL) - Alfa*CW2*MW*((q22-q32)*
     -        bup2*EL) + Alfa*MW*MW2*(C2B*bup1*EL)))/SW + Alfa/(32.D0*
     -        CW2*Pi*SW2)*(Bval(bb0,iintB10)*(MW*(C2B*bup1*EL) - 8*CW2*
     -        MW*(bup2*EL)))/SW + 1/(8.D0*CW2*Pi*SW2)*(Cval(cc1,iintB9)*
     -        (2*Alfa*CW2**2*MW*((q12-q22)*bup2*EL) - 2*Alfa*CW2*MW*SW2*
     -        (q32*bup2*EL) + 2*Alfa*CW2*MW*SW2*(q32*bup2*EL)))/SW 
         
         loopa1 = loopa1 + Cval(cc0,iintB8)*((5*Alfa*MW)/(8.D0*Pi*SW2)*
     -        ((q12+q22-q32)*bup2*EL)/SW - (Alfa*MW*MW2*(-2*CW2**2 + 
     -        SW2**2))/(4.D0*CW2**2*Pi*SW2)*(bup2*EL)/SW) + Alfa/(8.D0*
     -        CW2*Pi*SW2)*(Cval(cc00,iintB8)*(-(MW*(C2B*bup1*EL)) + 
     -        18*CW2*MW*(bup2*EL) + 2*MW*SW2*(bup2*EL)))/SW + Alfa/
     -        (8.D0*CW2*Pi*SW2)*(Cval(cc00,iintB9)*(-(MW*(C2B*bup1*EL))+ 
     -        CW2*(18*CW2*MW*(bup2*EL) + 2*MW*SW2*(bup2*EL))))/SW + 
     -        1/(8.D0*CW2**2*Pi*SW2)*(Cval(cc0,iintB9)*(4*Alfa*CW2**3*
     -        MW*((q12+1/2.D0*q22-q32)*bup2*EL) + 2*Alfa*CW2**2*MW*SW2*
     -        ((q12-1/2.D0*q22-q32)*bup2*EL) - Alfa*CW2**2*MW*SW2*((q12-
     -        2*q22)*bup2*EL) + Alfa*(MW*MW2*SW2**2*(C2B*bup1*EL) + 
     -        CW2**2*(4*CW2*MW*MZ2*(bup2*EL) - 2*MW*MZ2*SW2*(bup2*
     -        EL)))))/SW - (Alfa*MW)/(16.D0*CW2*Pi*SW2)*(EL*
     -        Bval(bb0,iintB6)*(2*CW2*(SBA*ZHggs(HiggsType,1)) + 2*CW2*
     -        (CBA*ZHggs(HiggsType,2)) + C2B*(SAB*ZHggs(HiggsType,1) - 
     -        CAB*ZHggs(HiggsType,2))))/SW


         loopa2 = loopa2-(2*Alfa*MW)/Pi*(bup2*EL*Cval(cc12,iintB5))/SW + 
     -        Alfa/(8.D0*CW2*Pi*SW2)*(Cval(cc12,iintB8)*(MW*(C2B*bup1*
     -        EL) - 18*CW2*MW*(bup2*EL) - 2*MW*SW2*(bup2*EL)))/SW + 
     -        1/(8.D0*CW2*Pi*SW2)*(Cval(cc1,iintB8)*(-(2*Alfa*MW*SW2*
     -        (bup2*EL)) + Alfa*CW2*MW*((-2)*bup2*EL)))/SW + 
     -        1/(8.D0*CW2*Pi*SW2)*(Cval(cc2,iintB8)*(-(2*Alfa*MW*SW2*
     -        (bup2*EL)) + 3*Alfa*CW2*MW*((-2/3D0)*bup2*EL)))/SW + 
     -        1/(8.D0*Pi)*(Cval(cc1,iintB5)*(Alfa*MW*((-2)*bup2*EL) + 
     -        Alfa*MW*((-2)*bup2*EL) - Alfa*MW*(4*bup2*EL)))/SW + 
     -        1/(8.D0*Pi*SW2)*(Cval(cc2,iintB9)*(-(Alfa*CW2*MW*((-2)*
     -        bup2*EL)) + Alfa*MW*SW2*((-2)*bup2*EL) - Alfa*MW*SW2*(4*
     -        bup2*EL)))/SW + 1/(8.D0*Pi)*(Cval(cc2,iintB5)*(3*Alfa*MW*
     -        ((-2/3D0)*bup2*EL) + Alfa*MW*((-4)*bup2*EL) - Alfa*MW*(2*
     -        bup2*EL)))/SW + 1/(8.D0*CW2*Pi)*(Cval(cc0,iintB5)*
     -        (5*Alfa*CW2*MW*((8/5D0)*bup2*EL)-Alfa*CW2*MW*(4*bup2*EL) - 
     -        Alfa*CW2*MW*(4*bup2*EL)))/SW + 1/(8.D0*CW2*Pi*SW2)*
     -        (Cval(cc1,iintB9)*(2*Alfa*CW2**2*MW*(10*bup2*EL) - 
     -        2*Alfa*CW2*MW*SW2*((-3/2D0)*bup2*EL) + 
     -        2*Alfa*CW2*MW*SW2*((-1/2D0)*bup2*EL) + 
     -        Alfa*((-(MW*(C2B*bup1*EL)) -2*CW2**2*MW*(bup2*EL)))))/SW + 
     -        Alfa/(8.D0*CW2*Pi*SW2)*((Cval(cc11,iintB9) + 
     -        Cval(cc12,iintB9))*(-(MW*(C2B*bup1*EL)) + 
     -        CW2*(18*CW2*MW*(bup2*EL) + 2*MW*SW2*(bup2*EL))))/SW
         
         loopa2 = loopa2 + Cval(cc0,iintB8)*((5*Alfa*MW)/(8.D0*Pi*SW2)*
     -        ((8/5D0)*bup2*EL)/SW) + 1/(8.D0*CW2**2*Pi*SW2)*
     -        (Cval(cc0,iintB9)*(4*Alfa*CW2**3*MW*((5/2D0)*bup2*EL) + 
     -        2*Alfa*CW2**2*MW*SW2*(bup2*EL)))/SW 


** neutral higgs
         do Hig4 = 1,3
            bup3 = SAB*UHggs(Hig4,1) - CAB*UHggs(Hig4,2)
            bup4 = SBA*UHggs(Hig4,1) + CBA*UHggs(Hig4,2)
            bup5 = -(CBA*UHggs(Hig4,1)) + SBA*UHggs(Hig4,2)
            bup6 = cI*(CBA*UHggs(Hig4,1)) - 
     -           cI*(SBA*UHggs(Hig4,2)) + UHggs(Hig4,3)
            bup7 = -(cI*(CBA*UHggs(Hig4,1))) + 
     -           cI*(SBA*UHggs(Hig4,2)) + UHggs(Hig4,3)
            bup8 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
            bup9 = SBA*ZHggs(HiggsType,1) + CBA*ZHggs(HiggsType,2)
            bup10 = 2*CW2*(SBA*ZHggs(HiggsType,1)) + 
     -           2*CW2*(CBA*ZHggs(HiggsType,2)) + 
     -           C2B*(SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2))
            bup11 = 2*CW2*(SBA*ZHggs(HiggsType,1)) + 
     -           2*CW2*(CBA*ZHggs(HiggsType,2)) + C2B*(-(SAB*
     -           ZHggs(HiggsType,1)) + CAB*ZHggs(HiggsType,2))
            bup12 = CBA*ZHggs(HiggsType,1) - SBA*ZHggs(HiggsType,2)

            loopa1 = loopa1 - (Alfa*cI*MW)/(8.D0*CW2*Pi*SW2)*(bup4*bup6*
     -           bup8*EL*S2B*Cval(cc00,iintB11(Hig4)))/SW+(Alfa*cI*MW)/
     -           (8.D0*CW2*Pi*SW2)*(bup4*bup7*bup8*EL*S2B*
     -           Cval(cc00,iintB12(Hig4)))/SW+(Alfa*MW)/(8.D0*CW2*Pi*
     -           SW2)*(bup11*bup4**2*EL*Cval(cc00,iintB13(Hig4)))/SW + 
     -           (Alfa*MW)/(8.D0*CW2*Pi*SW2)*(bup10*bup6*bup7*EL*
     -           Cval(cc00,iintB14(Hig4)))/SW + (Alfa*MW)/(4.D0*Pi*SW2)*
     -           (bup4**2*bup9*EL*Cval(cc0,iintB13(Hig4))*(-MW2))/SW - 
     -           (Alfa*cI*MW)/(8.D0*Pi*SW2)*((q12-q22-q32)*bup4*EL*
     -           Cval(cc2,iintB15(Hig4))*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))/SW - (Alfa*cI*MW)/(8.D0*Pi*SW2)*
     -           ((q12-q22+q32)*bup4*EL*Cval(cc2,iintB16(Hig4))*
     -           (bup12*UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3)))/SW - 
     -           (Alfa*MW)/(8.D0*Pi*SW2)*(bup4*EL*
     -           (Bval(bb0,iintB17(Hig4)) + Bval(bb0,iintB18(Hig4)))*
     -           (UHggs(Hig4,1)*ZHggs(HiggsType,1) + UHggs(Hig4,2)*
     -           ZHggs(HiggsType,2) + UHggs(Hig4,3)*
     -           ZHggs(HiggsType,3)))/SW
            
            loopa1 = loopa1 + Cval(cc0,iintB16(Hig4))*((Alfa*cI*MW)/
     -           (8.D0*Pi*SW2)*((q12-q32)*bup4*EL*(bup12*UHggs(Hig4,3)+ 
     -           bup5*ZHggs(HiggsType,3)))/SW - (Alfa*cI*MW*MW2)/(8.D0*
     -           Pi*SW2)*(bup4*EL*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))/SW) + Cval(cc0,iintB15(Hig4))*
     -           (-((Alfa*cI*MW)/(8.D0*Pi*SW2)*((q22-q32)*bup4*EL*
     -           (bup12*UHggs(Hig4,3)+ bup5*ZHggs(HiggsType,3)))/SW) + 
     -           (Alfa*cI*MW*MW2)/(8.D0*Pi*SW2)*(bup4*EL*(bup12*
     -           UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3)))/SW) + 
     -           Cval(cc1,iintB16(Hig4))*(-((Alfa*cI*MW)/(8.D0*Pi*SW2)*
     -           ((q12-q22-q32)*bup4*EL*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))/SW)) + Cval(cc1,iintB15(Hig4))*
     -           (-((Alfa*cI*MW)/(8.D0*Pi*SW2)*((q12-q22+q32)*bup4*EL*
     -           (bup12*UHggs(Hig4,3)+bup5*ZHggs(HiggsType,3)))/SW)) + 
     -           Cval(cc00,iintB15(Hig4))*(-((Alfa*cI*MW)/(8.D0*CW2*Pi*
     -           SW2)*(bup4*EL*S2B*(bup8*UHggs(Hig4,3) + bup3*
     -           ZHggs(HiggsType,3)))/SW) - (Alfa*cI*MW)/(8.D0*Pi*SW2)*
     -           (bup4*EL*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))/SW - (Alfa*cI*MW)/(8.D0*CW2*Pi*
     -           SW2)*(bup4*EL*SW*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3))))

            loopa1 = loopa1 + Cval(cc00,iintB16(Hig4))*
     -           ((Alfa*cI*MW)/(8.D0*CW2*Pi*SW2)*(bup4*EL*S2B*
     -           (bup8*UHggs(Hig4,3) + bup3*ZHggs(HiggsType,3)))/SW + 
     -           (Alfa*cI*MW)/(8.D0*Pi*SW2)*(bup4*EL*
     -           (bup12*UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3)))/SW + 
     -           (Alfa*cI*MW)/(8.D0*CW2*Pi*SW2)*(bup4*EL*SW*
     -           (bup12*UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3))))
            
            
            loopa2 = loopa2 + (Alfa*MW)/(8.D0*CW2*Pi*SW2)*(bup11*
     -           bup4**2*EL*Cval(cc1,iintB13(Hig4)))/SW - (Alfa*cI*MW)/
     -           (8.D0*CW2*Pi*SW2)*(bup4*bup6*bup8*EL*S2B*
     -           (Cval(cc11,iintB11(Hig4)) + Cval(cc12,iintB11(Hig4))))/
     -           SW + (Alfa*cI*MW)/(8.D0*CW2*Pi*SW2)*(bup4*bup7*bup8*EL*
     -           S2B*(Cval(cc11,iintB12(Hig4)) + 
     -           Cval(cc12,iintB12(Hig4))))/SW + (Alfa*MW)/(8.D0*CW2*Pi*
     -           SW2)*(bup11*bup4**2*EL*(Cval(cc11,iintB13(Hig4)) + 
     -           Cval(cc12,iintB13(Hig4))))/SW + (Alfa*MW)/(8.D0*CW2*Pi*
     -           SW2)*(bup10*bup6*bup7*EL*(Cval(cc1,iintB14(Hig4)) + 
     -           Cval(cc11,iintB14(Hig4)) + Cval(cc12,iintB14(Hig4))))/
     -           SW + (Alfa*MW)/(4.D0*Pi*SW2)*(bup4**2*bup9*EL*
     -           Cval(cc2,iintB13(Hig4)))/SW + (Alfa*MW)/(4.D0*Pi*SW2)*
     -           (bup4**2*bup9*EL*Cval(cc0,iintB13(Hig4)))/SW +(Alfa*cI*
     -           MW)/(8.D0*CW2*Pi*SW2)*(bup4*bup7*EL*
     -           Cval(cc1,iintB12(Hig4))*(-(2*CW2*(CBA*
     -           ZHggs(HiggsType,1))) + 2*CW2*(SBA*
     -           ZHggs(HiggsType,2)) + 2*cI*CW2*ZHggs(HiggsType,3) +
     -           S2B*SAB*ZHggs(HiggsType,1) - CAB*S2B*
     -           ZHggs(HiggsType,2)))/SW + (Alfa*cI*MW)/(8.D0*CW2*Pi*
     -           SW2)*(bup4*bup6*EL*Cval(cc1,iintB11(Hig4))*(2*CW2*(CBA*
     -           ZHggs(HiggsType,1)) -2*CW2*(SBA*ZHggs(HiggsType,2)) + 
     -           2*cI*CW2*ZHggs(HiggsType,3) - S2B*SAB*
     -           ZHggs(HiggsType,1) +CAB*S2B*ZHggs(HiggsType,2)))/SW - 
     -           (Alfa*cI*MW)/(8.D0*Pi*SW2)*(4*bup4*EL*
     -           Cval(cc2,iintB15(Hig4))*(bup12*UHggs(Hig4,3)+ bup5*
     -           ZHggs(HiggsType,3)))/SW-(Alfa*cI*MW)/(8.D0*Pi*SW2)*
     -           ((-4)*bup4*EL*Cval(cc2,iintB16(Hig4))*
     -           (bup12*UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3)))/SW 

            loopa2 = loopa2 + Cval(cc0,iintB16(Hig4))*((Alfa*cI*MW)/
     -           (8.D0*Pi*SW2)*(4*bup4*EL*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))/SW) + Cval(cc0,iintB15(Hig4))*
     -           (-((Alfa*cI*MW)/(8.D0*Pi*SW2)*(4*bup4*EL*
     -           (bup12*UHggs(Hig4,3)+bup5*ZHggs(HiggsType,3)))/SW)) + 
     -           Cval(cc12,iintB16(Hig4))*(-((Alfa*cI*MW)/(8.D0*CW2*Pi*
     -           SW2)*(bup4*EL*S2B*(bup8*UHggs(Hig4,3) + bup3*
     -           ZHggs(HiggsType,3)))/SW) - (Alfa*cI*MW)/(8.D0*Pi*SW2)*
     -           (bup4*EL*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))/SW - (Alfa*cI*MW)/(8.D0*CW2*Pi*SW2)*
     -           (bup4*EL*SW*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))) + Cval(cc12,iintB15(Hig4))*
     -           ((Alfa*cI*MW)/(8.D0*CW2*Pi*SW2)*(bup4*EL*S2B*
     -           (bup8*UHggs(Hig4,3) + bup3*ZHggs(HiggsType,3)))/SW + 
     -           (Alfa*cI*MW)/(8.D0*Pi*SW2)*(bup4*EL*
     -           (bup12*UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3)))/SW + 
     -           (Alfa*cI*MW)/(8.D0*CW2*Pi*SW2)*(bup4*EL*SW*
     -           (bup12*UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3)))) + 
     -           Cval(cc1,iintB16(Hig4))*(-((Alfa*cI*MW)/(8.D0*Pi*SW2)*
     -           ((-2)*bup4*EL*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))/SW) -(Alfa*cI*MW)/(4.D0*CW2*Pi*
     -           SW2)*(bup4*EL*SW*(bup12*UHggs(Hig4,3) + bup5*
     -           ZHggs(HiggsType,3)))) + Cval(cc1,iintB15(Hig4))*
     -           (-((Alfa*cI*MW)/(8.D0*Pi*SW2)*(2*bup4*EL*
     -           (bup12*UHggs(Hig4,3) +bup5*ZHggs(HiggsType,3)))/SW) + 
     -           (Alfa*cI*MW)/(4.D0*CW2*Pi*SW2)*(bup4*EL*SW*
     -           (bup12*UHggs(Hig4,3) + bup5*ZHggs(HiggsType,3)))) 
         enddo


         do Hig5 = 1,3
            do Hig4 = 1,3
               bup13 = SBA*UHggs(Hig4,1) + CBA*UHggs(Hig4,2)
               bup14 = -(cI*UHggs(Hig4,3)) + CBA*UHggs(Hig4,1) - 
     -              SBA*UHggs(Hig4,2)
               bup15 = SAB*UHggs(Hig5,1) - CAB*UHggs(Hig5,2)
               bup16 = SBA*UHggs(Hig5,1) + CBA*UHggs(Hig5,2)
               bup17 = cI*UHggs(Hig5,3) + CBA*UHggs(Hig5,1) - 
     -              SBA*UHggs(Hig5,2)
               bup18 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
               bup19 = 2*CAB*S2A*ZHggs(HiggsType,1) + C2A*(SAB*
     -              ZHggs(HiggsType,1) - 3*CAB*ZHggs(HiggsType,2))
               bup20 = -2*S2A*SAB*ZHggs(HiggsType,1) + 2*CAB*S2A*
     -              ZHggs(HiggsType,2) + C2A*(CAB*ZHggs(HiggsType,1) + 
     -              SAB*ZHggs(HiggsType,2))
               bup21 = 2*S2A*SAB*ZHggs(HiggsType,2) + C2A*(3*SAB*
     -              ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2))
               bup22 = C2B*UHggs(Hig4,3)*
     -              (bup18*UHggs(Hig5,3) + bup15*ZHggs(HiggsType,3)) - 
     -              UHggs(Hig4,2)*(bup20*UHggs(Hig5,1) + 
     -              bup19*UHggs(Hig5,2) + 
     -              C2B*CAB*UHggs(Hig5,3)*ZHggs(HiggsType,3)) + 
     -              UHggs(Hig4,1)*(bup21*UHggs(Hig5,1) - 
     -              bup20*UHggs(Hig5,2) + 
     -              C2B*SAB*UHggs(Hig5,3)*ZHggs(HiggsType,3))

               loopa1 = loopa1 - (Alfa*MW*MW2)/(8.D0*CW2*Pi*SW2)*
     -              (bup13*bup16*bup22*EL*Cval(cc0,iintB19(Hig4,Hig5)))/
     -              SW + (Alfa*MW)/(8.D0*CW2*Pi*SW2)*(bup13*bup16*bup22*
     -              EL*Cval(cc00,iintB19(Hig4,Hig5)))/SW + (Alfa*MW)/
     -              (8.D0*CW2*Pi*SW2)*(bup14*bup17*bup22*EL*
     -              Cval(cc00,iintB20(Hig4,Hig5)))/SW - (Alfa*MW)/
     -              (32.D0*CW2*Pi*SW2)*(bup22*EL*
     -              Bval(bb0,iintB21(Hig4,Hig5))*(UHggs(Hig4,1)*
     -              UHggs(Hig5,1) + UHggs(Hig4,2)*UHggs(Hig5,2) + 
     -              UHggs(Hig4,3)*UHggs(Hig5,3)))/SW


               loopa2 = loopa2 - (Alfa*MW)/(8.D0*CW2*Pi*SW2)*(bup13*
     -              bup16*bup22*EL*Cval(cc12,iintB19(Hig4,Hig5)))/
     -              SW - (Alfa*MW)/(8.D0*CW2*Pi*SW2)*(bup14*bup17*bup22*
     -              EL*Cval(cc12,iintB20(Hig4,Hig5)))/SW 
            enddo
         enddo


** chargino and neutralinos
         if (sector .ge. 4) then
            do Neu4 = 1,4
               do Cha5 = 1,2
                  do Cha4 = 1,2
                     bup37 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -                    CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                     bup38 = cI*(CB*ZHggs(HiggsType,3)) + 
     -                    CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                     bup39 = -(cI*(CA*ZHggs(HiggsType,1))) - 
     -                    cI*(SA*ZHggs(HiggsType,2)) + CB*ZHggs(HiggsType,3)
                     bup40 = cI*(CA*ZHggs(HiggsType,1)) + 
     -                    cI*(SA*ZHggs(HiggsType,2)) + 
     -                    CB*ZHggs(HiggsType,3)
                     bup41 = cI*(SB*ZHggs(HiggsType,3)) + 
     -                    SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                     bup42 = cI*(SB*ZHggs(HiggsType,3)) - 
     -                    SA*ZHggs(HiggsType,1) + CA*ZHggs(HiggsType,2)
                     bup43 = cI*(SA*ZHggs(HiggsType,1)) - cI*
     -                    (CA*ZHggs(HiggsType,2)) + 
     -                    SB*ZHggs(HiggsType,3)
                     bup44 = -(cI*(SA*ZHggs(HiggsType,1))) + 
     -                    cI*(CA*ZHggs(HiggsType,2)) + 
     -                    SB*ZHggs(HiggsType,3)
                     bup45 = 2*UChaC(Cha5,1)*ZNeu(Neu4,2) + 
     -                    Sqrt(2.D0)*UChaC(Cha5,2)*ZNeu(Neu4,3)
                     bup46 = -2*VChaC(Cha4,1)*ZNeu(Neu4,2) + 
     -                    Sqrt(2.D0)*VChaC(Cha4,2)*ZNeu(Neu4,4)
                     bup47 = 2*UCha(Cha4,1)*ZNeuC(Neu4,2) + 
     -                    Sqrt(2.D0)*UCha(Cha4,2)*ZNeuC(Neu4,3)
                     bup48 = -2*VCha(Cha5,1)*ZNeuC(Neu4,2) + 
     -                    Sqrt(2.D0)*VCha(Cha5,2)*ZNeuC(Neu4,4)
                     bup49 = cI*(bup45*bup47*(bup41*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup37*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2))) + bup46*bup48*(bup43*
     -                    UCha(Cha5,2)*VCha(Cha4,1) + bup39*
     -                    UCha(Cha5,1)*VCha(Cha4,2))
                     bup50 = cI*(bup46*bup48*(bup41*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup37*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2))) + bup45*bup47*(bup43*
     -                    UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup39*UCha(Cha5,1)*VCha(Cha4,2))
                     bup52 = cI*(bup45*bup46*
     -                    (bup42*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup38*UCha(Cha5,1)*VCha(Cha4,2))) + bup47*
     -                    bup48*(bup44*UChaC(Cha4,2)*VChaC(Cha5,1) + 
     -                    bup40*UChaC(Cha4,1)*VChaC(Cha5,2))

                     bup51 = cI*((q12+q22-q32)*bup46*bup48*
     -                    (bup41*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    bup37*UChaC(Cha4,1)*VChaC(Cha5,2))) + 
     -                    (q12+q22-q32)*bup45*bup47*
     -                    (bup43*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup39*UCha(Cha5,1)*VCha(Cha4,2))


                     loopa1 = loopa1 + (Alfa*cI)/(4.D0*Pi*SW2)*
     -                    (Cval(cc00,iintB24(Cha4,Cha5,Neu4))*
     -                    (-(bup50*EL*MCha(Cha4)) - bup49*EL*
     -                    MCha(Cha5)))/(SW*Sqrt(2.D0)) + (Alfa*cI)/
     -                    (8.D0*Pi*SW2)*(Bval(bb0,iintB25(Cha4,Cha5))*
     -                    (bup50*EL*MCha(Cha4) + bup49*EL*MCha(Cha5) - 
     -                    bup52*EL*MNeu(Neu4)))/(SW*Sqrt(2.D0)) + 
     -                    Cval(cc1,iintB24(Cha4,Cha5,Neu4))*((Alfa*cI)/
     -                    (8.D0*Pi*SW2)*((q12-q22)*bup52*EL*MNeu(Neu4))/
     -                    (SW*Sqrt(2.D0)) - (Alfa*cI)/(16.D0*Pi*SW2)*
     -                    (EL*MCha(Cha4)*(cI*((q12-q22-q32)*bup46*bup48*
     -                    (bup41*UChaC(Cha4,2)*VChaC(Cha5,1) - bup37*
     -                    UChaC(Cha4,1)*VChaC(Cha5,2))) + (q12-q22-q32)*
     -                    bup45*bup47*(bup43*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup39*UCha(Cha5,1)*VCha(Cha4,2))))/(SW*
     -                    Sqrt(2.D0)) - (Alfa*cI)/(16.D0*Pi*SW2)*(EL*
     -                    MCha(Cha5)*(cI*((q12-q22+q32)*bup45*bup47*
     -                    (bup41*UChaC(Cha4,2)*VChaC(Cha5,1) - bup37*
     -                    UChaC(Cha4,1)*VChaC(Cha5,2))) + (q12-q22+q32)*
     -                    bup46*bup48*(bup43*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup39*UCha(Cha5,1)*VCha(Cha4,2))))/(SW*
     -                    Sqrt(2.D0))) + 
     -                    Cval(cc0,iintB24(Cha4,Cha5,Neu4))*(-((Alfa*
     -                    cI)/(16.D0*Pi*SW2)*(bup51*EL*MCha(Cha4))/(SW*
     -                    Sqrt(2.D0))) - (Alfa*cI)/(8.D0*Pi*SW2)*(q12*
     -                    bup49*EL*MCha(Cha5))/(SW*Sqrt(2.D0)) + (Alfa*
     -                    cI)/(8.D0*Pi*SW2)*(q12*bup52*EL*MNeu(Neu4))/
     -                    (SW*Sqrt(2.D0)) + (Alfa*cI)/(8.D0*Pi*SW2)*
     -                    (bup50*EL*MCha(Cha4)*MNeu2(Neu4))/(SW*
     -                    Sqrt(2.D0))+(Alfa*cI)/(8.D0*Pi*SW2)*(bup49*EL*
     -                    MCha(Cha5)*MNeu2(Neu4))/(SW*Sqrt(2.D0))-
     -                    (Alfa*cI)/(8.D0*Pi*SW2)*(bup52*EL*MNeu(Neu4)*
     -                    MNeu2(Neu4))/(SW*Sqrt(2.D0))-(Alfa*cI)/(8.D0*
     -                    Pi*SW2)*(EL*MCha(Cha4)*MCha(Cha5)*MNeu(Neu4)*
     -                    (cI*(bup47*bup48*(bup42*UCha(Cha5,2)*
     -                    VCha(Cha4,1) + bup38*UCha(Cha5,1)*
     -                    VCha(Cha4,2))) + bup45*bup46*(bup44*
     -                    UChaC(Cha4,2)*VChaC(Cha5,1) + bup40*
     -                    UChaC(Cha4,1)*VChaC(Cha5,2))))/(SW*
     -                    Sqrt(2.D0)))
                     
                     loopa1 = loopa1 +Cval(cc2,iintB24(Cha4,Cha5,Neu4))*
     -                    (-((Alfa*cI)/(16.D0*Pi*SW2)*(bup51*EL*
     -                    MCha(Cha4))/(SW*Sqrt(2.D0))) -(Alfa*cI)/(8.D0*
     -                    Pi*SW2)*(q12*bup49*EL*MCha(Cha5))/(SW*
     -                    Sqrt(2.D0)) + (3*Alfa*cI)/(16.D0*Pi*SW2)*(EL*
     -                    MNeu(Neu4)*(cI*((q12+1/3.D0*q22-1/3.D0*q32)*
     -                    bup45*bup46*(bup42*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup38*UCha(Cha5,1)*VCha(Cha4,2))) + 
     -                    (q12+1/3.D0*q22-1/3.D0*q32)*bup47*bup48*
     -                    (bup44*UChaC(Cha4,2)*VChaC(Cha5,1) + 
     -                    bup40*UChaC(Cha4,1)*VChaC(Cha5,2))))/
     -                    (SW*Sqrt(2.D0)))


                     bup51 = cI*(2*bup46*bup48*(bup41*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup37*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2))) + 2*bup45*bup47*(bup43*
     -                    UCha(Cha5,2)*VCha(Cha4,1) +bup39*UCha(Cha5,1)*
     -                    VCha(Cha4,2))
                     
                     loopa2 = loopa2 +Cval(cc2,iintB24(Cha4,Cha5,Neu4))*
     -                    (-((Alfa*cI)/(16.D0*Pi*SW2)*(bup51*EL*
     -                    MCha(Cha4))/(SW*Sqrt(2.D0))) + (3*Alfa*cI)/
     -                    (16.D0*Pi*SW2)*(EL*MNeu(Neu4)*(cI*((2/3D0)*
     -                    bup45*bup46*(bup42*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup38*UCha(Cha5,1)*VCha(Cha4,2))) + (2/3D0)*
     -                    bup47*bup48*(bup44*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1)+ bup40*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2))))/(SW*Sqrt(2.D0))) + 
     -                    Cval(cc1,iintB24(Cha4,Cha5,Neu4))*(-(Alfa*cI)/
     -                    (16.D0*Pi*SW2)*(EL*MCha(Cha4)*(cI*(6*bup46*
     -                    bup48*(bup41*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    bup37*UChaC(Cha4,1)*VChaC(Cha5,2))) + 6*bup45*
     -                    bup47*(bup43*UCha(Cha5,2)*VCha(Cha4,1) +bup39*
     -                    UCha(Cha5,1)*VCha(Cha4,2))))/(SW*Sqrt(2.D0)) - 
     -                    (Alfa*cI)/(16.D0*Pi*SW2)*(EL*MCha(Cha5)*(cI*
     -                    (2*bup45*bup47*(bup41*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup37*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2))) + 2*bup46*bup48*(bup43*
     -                    UCha(Cha5,2)*VCha(Cha4,1) +bup39*UCha(Cha5,1)*
     -                    VCha(Cha4,2))))/(SW*Sqrt(2.D0))) + (Alfa*cI)/
     -                    (4.D0*Pi*SW2)*
     -                    ((Cval(cc11,iintB24(Cha4,Cha5,Neu4)) + 
     -                    Cval(cc12,iintB24(Cha4,Cha5,Neu4)))*(-(bup50*
     -                    EL*MCha(Cha4)) - bup49*EL*MCha(Cha5)))/(SW*
     -                    Sqrt(2.D0))+Cval(cc0,iintB24(Cha4,Cha5,Neu4))*
     -                    (-((Alfa*cI)/(16.D0*Pi*SW2)*(bup51*EL*
     -                    MCha(Cha4))/(SW*Sqrt(2.D0)))) 


                     bup51 = cI*(2*bup46*bup48*(bup41*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup37*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2))) + (-2)*bup45*bup47*
     -                    (bup43*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup39*UCha(Cha5,1)*VCha(Cha4,2))

                     loopa3 = loopa3 +Cval(cc1,iintB24(Cha4,Cha5,Neu4))*
     -                    (-(Alfa*cI)/(16.D0*Pi*SW2)*(EL*MCha(Cha4)*
     -                    (cI*(2*bup46*bup48*(bup41*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup37*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2))) + (-2)*bup45*bup47*
     -                    (bup43*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup39*UCha(Cha5,1)*VCha(Cha4,2))))/
     -                    (SW*Sqrt(2.D0)) - (Alfa*cI)/(16.D0*Pi*SW2)*
     -                    (EL*MCha(Cha5)*(cI*((-2)*bup45*bup47*
     -                    (bup41*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    bup37*UChaC(Cha4,1)*VChaC(Cha5,2))) + 2*bup46*
     -                    bup48*(bup43*UCha(Cha5,2)*VCha(Cha4,1)+ bup39*
     -                    UCha(Cha5,1)*VCha(Cha4,2))))/(SW*Sqrt(2.D0)))+ 
     -                    Cval(cc0,iintB24(Cha4,Cha5,Neu4))*
     -                    (-((Alfa*cI)/(16.D0*Pi*SW2)*
     -                    (bup51*EL*MCha(Cha4))/(SW*Sqrt(2.D0))))
                     
                     loopa3 = loopa3 +Cval(cc2,iintB24(Cha4,Cha5,Neu4))*
     -                    (-((Alfa*cI)/(16.D0*Pi*SW2)*
     -                    (bup51*EL*MCha(Cha4))/(SW*Sqrt(2.D0))) + 
     -                    (3*Alfa*cI)/(16.D0*Pi*SW2)*(EL*MNeu(Neu4)*
     -                    (cI*((-2/3D0)*bup45*bup46*
     -                    (bup42*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup38*UCha(Cha5,1)*VCha(Cha4,2))) + 
     -                    (2/3D0)*bup47*bup48*
     -                    (bup44*UChaC(Cha4,2)*VChaC(Cha5,1) + 
     -                    bup40*UChaC(Cha4,1)*VChaC(Cha5,2))))/
     -                    (SW*Sqrt(2.D0)))
                  enddo
               enddo
            enddo
            
            do Neu5 = 1,4
               do Neu4 = 1,4
                  do Cha4 = 1,2
                     bup23 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -                    CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                     bup24 = cI*(CB*ZHggs(HiggsType,3)) + 
     -                    CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                     bup25 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -                    SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                     bup26 = cI*(SB*ZHggs(HiggsType,3)) + 
     -                    SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                     bup27 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                     bup28 = 2*VChaC(Cha4,1)*ZNeu(Neu4,2) - 
     -                    Sqrt(2.D0)*VChaC(Cha4,2)*ZNeu(Neu4,4)
                     bup29 = -(CW*ZNeu(Neu5,2)) + SW*ZNeu(Neu5,1)
                     bup30 = 2*UChaC(Cha4,1)*ZNeu(Neu5,2) + 
     -                    Sqrt(2.D0)*UChaC(Cha4,2)*ZNeu(Neu5,3)
                     bup31 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                     bup32 = 2*UCha(Cha4,1)*ZNeuC(Neu4,2) + 
     -                    Sqrt(2.D0)*UCha(Cha4,2)*ZNeuC(Neu4,3)
                     bup33 = -(CW*ZNeuC(Neu5,2)) + SW*ZNeuC(Neu5,1)
                     bup34 = -2*VCha(Cha4,1)*ZNeuC(Neu5,2) + 
     -                    Sqrt(2.D0)*VCha(Cha4,2)*ZNeuC(Neu5,4)
                     bup35 = bup25*bup29*ZNeu(Neu4,3) + 
     -                    bup24*bup29*ZNeu(Neu4,4) + bup27*
     -                    (bup25*ZNeu(Neu5,3) + bup24*ZNeu(Neu5,4))
                     bup36 = bup26*bup33*ZNeuC(Neu4,3) + 
     -                    bup23*bup33*ZNeuC(Neu4,4) + bup31*
     -                    (bup26*ZNeuC(Neu5,3) + bup23*ZNeuC(Neu5,4))

                     loopa1 = loopa1 + 1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc2,iintB22(Cha4,Neu4,Neu5))*(Alfa*
     -                    ((q12+3*q22-q32)*bup32*bup34*bup35*EL*
     -                    MCha(Cha4)) - Alfa*((q12+3*q22-q32)*bup28*
     -                    bup30*bup36*EL*MCha(Cha4)) + Alfa*(q22*(-2*
     -                    bup28*bup34*bup35*EL*MNeu(Neu4) + 2*bup30*
     -                    bup32*bup36*EL*MNeu(Neu4))) + Alfa*((q12+q22-
     -                    q32)*bup30*bup32*bup35*EL*MNeu(Neu5)) - 
     -                    Alfa*((q12+q22-q32)*bup28*bup34*bup36*EL*
     -                    MNeu(Neu5))))/SW + 1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc1,iintB22(Cha4,Neu4,Neu5))*
     -                    (3*Alfa*((q12+1/3.D0*q22-1/3.D0*q32)*bup32*
     -                    bup34*bup35*EL*MCha(Cha4)) - 3*Alfa*((q12+
     -                    1/3.D0*q22-1/3.D0*q32)*bup28*bup30*bup36*EL*
     -                    MCha(Cha4)) - Alfa*((q12+q22-q32)*bup28*bup34*
     -                    bup35*EL*MNeu(Neu4)) + Alfa*((q12+q22-q32)*
     -                    bup30*bup32*bup36*EL*MNeu(Neu4)) + 
     -                    Alfa*(q12*(2*bup30*bup32*bup35*EL*MNeu(Neu5) - 
     -                    2*bup28*bup34*bup36*EL*MNeu(Neu5)))))/SW + 
     -                    Alfa/(8.D0*CW*Pi*SW2)*
     -                    (Cval(cc00,iintB22(Cha4,Neu4,Neu5))*(bup30*
     -                    bup32*(-(bup36*EL*MNeu(Neu4)) - bup35*EL*
     -                    MNeu(Neu5))+bup28*bup34*(bup35*EL*MNeu(Neu4) + 
     -                    bup36*EL*MNeu(Neu5))))/SW + Alfa/(16.D0*CW*Pi*
     -                    SW2)*(Bval(bb0,iintB23(Neu4,Neu5))*(bup32*
     -                    bup34*bup35*EL*MCha(Cha4) +bup30*bup32*(bup36*
     -                    EL*MNeu(Neu4) + bup35*EL*MNeu(Neu5)) + bup28*
     -                    (-(bup30*bup36*EL*MCha(Cha4)) +bup34*(-(bup35*
     -                    EL*MNeu(Neu4)) - bup36*EL*MNeu(Neu5)))))/SW
                     
                     loopa1 = loopa1 + 1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc0,iintB22(Cha4,Neu4,Neu5))*
     -                    (Alfa*((q12+q22-q32)*bup32*bup34*bup35*EL*
     -                    MCha(Cha4)) - Alfa*((q12+q22-q32)*bup28*bup30*
     -                    bup36*EL*MCha(Cha4)) + Alfa*(bup30*bup32*
     -                    (2*bup36*EL*MCha2(Cha4)*MNeu(Neu4) + 2*bup35*
     -                    EL*MCha2(Cha4)*MNeu(Neu5)) + bup32*bup34*
     -                    (2*bup35*EL*MCha(Cha4)*MCha2(Cha4) + 
     -                    2*bup36*EL*MCha(Cha4)*MNeu(Neu4)*MNeu(Neu5)) + 
     -                    bup28*(bup34*(-2*bup35*EL*MCha2(Cha4)*
     -                    MNeu(Neu4) - 2*bup36*EL*MCha2(Cha4)*
     -                    MNeu(Neu5)) + bup30*(-2*bup36*EL*MCha(Cha4)*
     -                    MCha2(Cha4) -2*bup35*EL*MCha(Cha4)*MNeu(Neu4)*
     -                    MNeu(Neu5))))))/SW


                     loopa2 = loopa2 + 1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc2,iintB22(Cha4,Neu4,Neu5))*
     -                    (Alfa*(2*bup32*bup34*bup35*EL*MCha(Cha4)) - 
     -                    Alfa*(2*bup28*bup30*bup36*EL*MCha(Cha4)) + 
     -                    Alfa*(2*bup30*bup32*bup35*EL*MNeu(Neu5)) - 
     -                    Alfa*(2*bup28*bup34*bup36*EL*MNeu(Neu5))))/SW+
     -                    1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc1,iintB22(Cha4,Neu4,Neu5))*(3*Alfa*
     -                    ((2/3D0)*bup32*bup34*bup35*EL*MCha(Cha4)) - 
     -                    3*Alfa*((2/3D0)*bup28*bup30*bup36*EL*
     -                    MCha(Cha4)) - Alfa*(2*bup28*bup34*bup35*EL*
     -                    MNeu(Neu4)) + Alfa*(2*bup30*bup32*bup36*EL*
     -                    MNeu(Neu4))))/SW + Alfa/(8.D0*CW*Pi*SW2)*
     -                    (Cval(cc12,iintB22(Cha4,Neu4,Neu5))*
     -                    (bup30*bup32*(bup36*EL*MNeu(Neu4) + 
     -                    bup35*EL*MNeu(Neu5)) +bup28*bup34*(-(bup35*EL*
     -                    MNeu(Neu4)) - bup36*EL*MNeu(Neu5))))/SW 
                     
                     loopa2 = loopa2 + 1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc0,iintB22(Cha4,Neu4,Neu5))*
     -                    (Alfa*(2*bup32*bup34*bup35*EL*MCha(Cha4)) - 
     -                    Alfa*(2*bup28*bup30*bup36*EL*MCha(Cha4))))/SW


                     loopa3 = loopa3 + 1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc2,iintB22(Cha4,Neu4,Neu5))*
     -                    (Alfa*((-2)*bup32*bup34*bup35*EL*MCha(Cha4)) - 
     -                    Alfa*(2*bup28*bup30*bup36*EL*MCha(Cha4)) + 
     -                    Alfa*((-2)*bup30*bup32*bup35*EL*MNeu(Neu5)) - 
     -                    Alfa*(2*bup28*bup34*bup36*EL*MNeu(Neu5))))/SW+
     -                    1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc1,iintB22(Cha4,Neu4,Neu5))*(3*Alfa*
     -                    ((-2/3D0)*bup32*bup34*bup35*EL*MCha(Cha4)) -3*
     -                    Alfa*((2/3D0)*bup28*bup30*bup36*EL*MCha(Cha4)) - 
     -                    Alfa*((-2)*bup28*bup34*bup35*EL*MNeu(Neu4)) + 
     -                    Alfa*(2*bup30*bup32*bup36*EL*MNeu(Neu4))))/SW
                     
                     loopa3 = loopa3 + 1/(32.D0*CW*Pi*SW2)*
     -                    (Cval(cc0,iintB22(Cha4,Neu4,Neu5))*
     -                    (Alfa*((-2)*bup32*bup34*bup35*EL*MCha(Cha4)) - 
     -                    Alfa*(2*bup28*bup30*bup36*EL*MCha(Cha4))))/SW
                  enddo
               enddo
            enddo
         end if  ! end of neutralinos/charginos
         
      end if
** End of bosonic contributions 


** The eps part of the loop - converting between FC notation
** and standard notation (in FC Eps = -i*levi-civita)
      loopa3 = loopa3*(0.D0,-1.D0)


** Now we want the form factors themselves         
      a2vtxW = -loopa2
         
      a1vtxW = loopa1 + (loopa2*q1q2) + ct 
         
      a3vtxW = loopa3

      if (debug .ne. 0) then
         write(*,*)"  "
         write(*,*)"WWh vtx a1 =", a1vtxW
         write(*,*)"WWh vtx a2 =", a2vtxW
         write(*,*)"WWh vtx a3 =", a3vtxW
         write(*,*)"WWh FERMION contn to vtx a1 =", a1vtxWfm
         write(*,*)"WWh FERMION contn to vtx a2 =", a2vtxWfm
         write(*,*)"WWh FERMION contn to vtx a3 =", a3vtxWfm
      end if

      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the ZZ-higgs vertex
*** in SUSY

      subroutine zz_vtxSUSY(a1vtxZ,a2vtxZ,a3vtxZ,q12,q22,q32,q1q2,debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Dummy variables (q-> 'quark', l-> 'lepton', b -> 'boson'
      double complex qup1, qup2, qup3, qup4, qup5, qup6
      double complex qup7, qup8, qup9, qup10, qup11, qup12
      double complex qup13, qup14, qup15, qup16, qup17, qup18, qup19
      double complex qup20, qup21, qup22, qup23, qup24, qup25, qup26
      double complex qup27, qup28, qup29, qup30, qup31, qup32, qup33
      double complex qup34, qup35, qup36, qup37, qup38, qup39, qup40
      double complex lup1, lup2, lup3, lup4, lup5, lup6, lup7, lup8
      double complex lup9, lup10, lup11, lup12, lup13, lup14, lup15
      double complex lup16, lup17, lup18, lup19, lup20, lup21
      double complex bup1, bup2, bup3,  bup4, bup5, bup6, bup7 
      double complex bup8, bup9, bup10, bup11, bup12, bup13, bup14
      double complex bup15, bup16, bup17, bup18, bup19, bup20, bup21
      double complex bup22, bup23, bup24, bup25, bup26, bup27, bup28
      double complex bup29, bup30, bup31, bup32, bup33, bup34, bup35
      double complex bup36, bup37, bup38, bup39, bup40, bup41, bup42
      double complex bup43, bup44, bup45, bup46, bup47, bup48, bup49
      double complex bup50, bup51, bup52, bup53, bup54, bup55, bup56

** The arguments of the integrals (suffix Q='quark', L='lepton', B='boson')
      integer*8 iintQ1(3), iintQ2(3), iintQ3(3), iintQ4(3), iintQ5(3,2,2)
      integer*8 iintQ6(3,2,2), iintQ7(3,2,2,2), iintQ8(3,2,2,2)
      integer*8 iintL1(3), iintL2(3), iintL3(3), iintL4(3), iintL5(3,2,2)
      integer*8 iintL6(3,2,2,2)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6, iintB7
      integer*8 iintB8(4,4,4), iintB9(4,4), iintB10(3), iintB11(3)
      integer*8 iintB12(3), iintB13(3,3), iintB14(3,3), iintB15(3,3)
      integer*8 iintB16(3,3), iintB17(3,3,3), iintB18(2,2,2), iintB19(2,2)


** These are the contributions from the counter-term (g-mu-nu), the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, cth0, ctHH, loopa1, loopa2, loopa3


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita) 
** for ZZ-higgs.  
      double complex a1vtxZ, a2vtxZ, a3vtxZ


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions
      double complex ctFM, cth0FM, ctHHfm
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


** Do loop variables
      integer Gen4, Sfe4, Sfe5, Sfe6, g4ini
      integer Neu4, Neu5, Neu6, Cha4, Cha5, Cha6
      integer Hig4, Hig5, Hig6


** Higgs mixing matrices
      double complex UHggs(3,3), ZHggs(3,3)

** EM coupling, identity matrix
      double precision Alfa, Alfa2
      double precision IndexDelta(4,4)

      integer debug

** Model parameters
#include "mssm.inc"
#include "global.inc"      

** This is needed for looptools. 
#include "looptools.h"


** Set UHggs and ZHggs to 1 for debugging
      do Hig4 = 1, 3
         do Hig5 = 1, 3
            if (Hig4 .eq. Hig5) then
*               ZHggs(Hig4,Hig5) = (1D0,0D0)
               UHggs(Hig4,Hig5) = (1D0,0D0)
            else
*               ZHggs(Hig4,Hig5) = (0D0,0D0)
               UHggs(Hig4,Hig5) = (0D0,0D0)
            end if
         end do
      end do


      if (propLEVEL .ne. 2) then
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               if (Hig4 .eq. Hig5) then
                  ZHggs(Hig4,Hig5) = (1D0,0D0)
               else
                  ZHggs(Hig4,Hig5) = (0D0,0D0)
               end if
            end do
         end do
      else
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               ZHggs(Hig4,Hig5) = Hgs(Hig4,Hig5)
            end do
         end do
      end if


** Clearing the values of the formfactors
      a1vtxZ = (0D0,0D0)
      a2vtxZ = (0D0,0D0)
      a3vtxZ = (0D0,0D0)


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      do Neu4 = 1, 4
         do Neu5 = 1, 4
            if (Neu4 .eq. Neu5) then
               IndexDelta(Neu4,Neu5) = 1D0
            else
               IndexDelta(Neu4,Neu5) = 0D0
            end if
         end do
      end do


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions, 3 for SM-type particles
** and sfermions and 4, 5, 6 for all particles
      if (sector .eq. 1) then
         g4ini = 3
      else 
         g4ini = 1
      end if


** Now for the arguments of the integrals.  

** First the quarks and squarks
      do Gen4 = g4ini, 3
         iintQ1(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ2(Gen4) = Cget(q22,q32,q12,Mf2(4,Gen4),
     -        Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ3(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
         iintQ4(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),
     -        Mf2(3,Gen4),Mf2(3,Gen4))
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3
               iintQ5(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,4,Gen4),MSf2(Sfe5,4,Gen4))
               iintQ6(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,3,Gen4))
            enddo
         enddo
      enddo
      
      do Sfe6 = 1,2
         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = g4ini, 3  
                  iintQ7(Gen4,Sfe4,Sfe5,Sfe6) = Cget(q12,q32,q22,
     -                 MSf2(Sfe4,4,Gen4),MSf2(Sfe5,4,Gen4),
     -                 MSf2(Sfe6,4,Gen4))
                  iintQ8(Gen4,Sfe4,Sfe5,Sfe6) = Cget(q12,q32,q22,
     -                 MSf2(Sfe4,3,Gen4),MSf2(Sfe5,3,Gen4),
     -                 MSf2(Sfe6,3,Gen4))               
               enddo
            enddo
         enddo
      enddo

** The (s)leptonic integrals
      if (sector .ne. 1) then
         do Gen4 = 1,3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL2(Gen4) = Cget(q22,q32,q12,Mf2(2,Gen4),
     -           Mf2(2,Gen4),Mf2(2,Gen4))
            iintL3(Gen4) = Bget(q32,MSf2(1,1,Gen4),MSf2(1,1,Gen4))
            iintL4(Gen4) = Cget(q22,q32,q12,MSf2(1,1,Gen4),
     -           MSf2(1,1,Gen4),MSf2(1,1,Gen4))
         enddo

         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  iintL5(Gen4,Sfe4,Sfe5) = 
     -                 Bget(q32,MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4))
               enddo
            enddo
         enddo

         do Sfe6 = 1,2
            do Sfe5 = 1,2
               do Sfe4 = 1,2
                  do Gen4 = 1,3  
                     iintL6(Gen4,Sfe4,Sfe5,Sfe6) = Cget(q12,q32,q22,
     -                    MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4),
     -                    MSf2(Sfe6,2,Gen4))                
                  enddo
               enddo
            enddo
         enddo
      end if

** Bosonic integrals and their superpartners
      if (sector .ge. 3) then
         iintB1 = Bget(q32,MZ2,MZ2)
         iintB2 = Bget(q12,MW2,MW2)
         iintB3 = Bget(q22,MW2,MW2)
         iintB4 = Bget(q32,MHp2,MHp2)
         iintB5 = Cget(q22,q32,q12,MHp2,MHp2,MHp2)
         iintB6 = Bget(q32,MW2,MW2)
         iintB7 = Cget(q22,q32,q12,MW2,MW2,MW2)

** neutral higgs ..
         do Hig4 = 1,3
            iintB10(Hig4) = Cget(q32,q12,q22,MZ2,MZ2,MHiggs2(Hig4))
            iintB11(Hig4) = Bget(q12,MZ2,MHiggs2(Hig4))
            iintB12(Hig4) = Bget(q22,MZ2,MHiggs2(Hig4))
         enddo

         do Hig5 = 1,3
            do Hig4 = 1,3
               iintB13(Hig4,Hig5) = 
     -              Cget(q12,q22,q32,MZ2,MHiggs2(Hig4),MHiggs2(Hig5))
               iintB14(Hig4,Hig5) = 
     -              Cget(q22,q12,q32,MZ2,MHiggs2(Hig4),MHiggs2(Hig5))
               iintB15(Hig4,Hig5) = 
     -              Cget(q12,q32,q22,MZ2,MHiggs2(Hig4),MHiggs2(Hig5))
               iintB16(Hig4,Hig5) =Bget(q32,MHiggs2(Hig4),MHiggs2(Hig5))
            enddo
         enddo

         do Hig6 = 1,3
            do Hig5 = 1,3
               do Hig4 = 1,3
                  iintB17(Hig4,Hig5,Hig6) = 
     -                 Cget(q12,q32,q22,MHiggs2(Hig4),MHiggs2(Hig5),
     -                 MHiggs2(Hig6))
               enddo
            enddo
         enddo

** chargino/neutralino integrals
         if (sector .ge. 4) then
            do Cha6 = 1,2 
               do Cha5 = 1,2
                  iintB19(Cha5,Cha6) = Bget(q32,MCha2(Cha5),MCha2(Cha6))
               enddo
            enddo
            do Neu6 = 1,4
               do Neu5 = 1,4
                  iintB9(Neu5,Neu6) = Bget(q32,MNeu2(Neu5),MNeu2(Neu6))
               enddo
            enddo
            do Cha6 = 1,2
               do Cha5 = 1,2
                  do Cha4 = 1,2
                     iintB18(Cha4,Cha5,Cha6) = Cget(q12,q32,q22,
     -                    MCha2(Cha4),MCha2(Cha5),MCha2(Cha6))
                  enddo
               enddo
            enddo
            do Neu6 = 1,4
               do Neu5 = 1,4
                  do Neu4 = 1,4
                     iintB8(Neu4,Neu5,Neu6) = Cget(q12,q32,q22,
     -                    MNeu2(Neu4),MNeu2(Neu5),MNeu2(Neu6))
                  enddo
               enddo
            enddo
         end if
      end if


** Clearing the values of each contribution
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0) 
      ct = (0D0,0D0)   

** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
      cth0 = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(CW2*MW2*SW2*(CBA*dZh0HH) + 
     -     2*CW2*MW2*SW2*(CB*CBA*dTB1*SB) - 4*MW2*SW2*(dSW1*SBA*SW) + 
     -     CW2*(SW2*(dMWsq1*SBA) + 2*MW2*SW2*(dZet1*SBA) + 
     -     MW2*SW2*(dZh0*SBA) + 2*MW2*SW2*(dZZZ1*SBA) + 
     -     2*MW2*(dSW1*SBA*SW))))/SW
         
      ctHH = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(CW2*MW2*SW2*(dZh0HH*SBA) - 
     -     2*CW2*MW2*SW2*(CB*dTB1*SB*SBA) + CBA*(-(4*MW2*SW2*
     -     (dSW1*SW)) + CW2*(SW2*dMWsq1 + 2*MW2*SW2*dZet1 + 
     -     MW2*SW2*dZHH + 2*MW2*SW2*dZZZ1 + 2*MW2*(dSW1*SW)))))/SW  

      if (sector .eq. 4) then
         cth0 = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(CW2*MW2*SW2*(CBA*
     -        dZh0HHfull) + 2*CW2*MW2*SW2*(CB*CBA*dTB1full*SB) - 4*MW2*
     -        SW2*(dSW1full*SBA*SW) + CW2*(SW2*(dMWsq1full*SBA) + 2*MW2*
     -        SW2*(dZet1full*SBA) + MW2*SW2*(dZh0full*SBA) + 2*MW2*SW2*
     -        (dZZZ1full*SBA) + 2*MW2*(dSW1full*SBA*SW))))/SW
         
         ctHH = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(CW2*MW2*SW2*(dZh0HHfull*
     -        SBA) - 2*CW2*MW2*SW2*(CB*dTB1full*SB*SBA) + CBA*(-(4*MW2*
     -        SW2*(dSW1full*SW)) + CW2*(SW2*dMWsq1full + 2*MW2*SW2*
     -        dZet1full + MW2*SW2*dZHHfull + 2*MW2*SW2*dZZZ1full + 
     -        2*MW2*(dSW1full*SW)))))/SW  
      end if

      ct = cth0*ZHggs(HiggsType,1) + ctHH*ZHggs(HiggsType,2)

** fermion part of the counterterms
      cth0FM = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(CW2*MW2*SW2*(CBA*dZh0HHfm) + 
     -     2*CW2*MW2*SW2*(CB*CBA*dTB1fm*SB) - 4*MW2*SW2*(dSW1fm*SBA*SW) + 
     -     CW2*(SW2*(dMWsq1fm*SBA) + 2*MW2*SW2*(dZet1fm*SBA) + 
     -     MW2*SW2*(dZh0fm*SBA) + 2*MW2*SW2*(dZZZ1fm*SBA) + 
     -     2*MW2*(dSW1fm*SBA*SW))))/SW
         
      ctHHfm = MW/(2.D0*CW2**2*MW2*SW2)*(EL*(CW2*MW2*SW2*(dZh0HHfm*SBA) - 
     -     2*CW2*MW2*SW2*(CB*dTB1fm*SB*SBA) + CBA*(-(4*MW2*SW2*
     -     (dSW1fm*SW)) + CW2*(SW2*dMWsq1fm + 2*MW2*SW2*dZet1fm + 
     -     MW2*SW2*dZHHfm + 2*MW2*SW2*dZZZ1fm + 2*MW2*(dSW1fm*SW)))))/SW 

      ctFM = cth0fm*ZHggs(HiggsType,1) + ctHHfm*ZHggs(HiggsType,2)

      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'ZZh ct=', ct
         write(*,*)'ZZh FERMION contn to ct=', ctFM
      end if
         

** Now the loop

** Quark and squark contributions
      do Gen4 = g4ini, 3
         qup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
         qup2 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
         qup3 = -(cI*(SB*ZHggs(HiggsType,3))) + CA*TB*
     -        ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)
         qup4 = cI*(SB*ZHggs(HiggsType,3))+ CA*TB*ZHggs(HiggsType,1) + 
     -        SA*TB*ZHggs(HiggsType,2)

         qup5 = q12*(4*SA*ZHggs(HiggsType,1)-4*CA*ZHggs(HiggsType,2))+ 
     -        (q12+q22-q32)*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -        (q12+q22-q32)*(cI*(SB*ZHggs(HiggsType,3)) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))
         qup6 = q22*(4*SA*ZHggs(HiggsType,1)-4*CA*ZHggs(HiggsType,2))+ 
     -        (q12+q22-q32)*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -        (q12+q22-q32)*(cI*(SB*ZHggs(HiggsType,3)) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))

         loopa1 = loopa1 + 1/(32.D0*CW2*MW*Pi*SW2)*
     -        (Cval(cc2,iintQ4(Gen4))*
     -        (3*Alfa*((q12+q22-q32)*qup3*EL*Mf2(3,Gen4)) + 
     -        3*Alfa*((q12+q22-q32)*qup4*EL*Mf2(3,Gen4)) + 
     -        12*Alfa*(q12*qup2*EL*TB*Mf2(3,Gen4))))/(SB*SW*TB)+
     -        1/(32.D0*CW2*MW*Pi*SW2)*(Cval(cc1,iintQ4(Gen4))*
     -        (3*Alfa*((q12+q22-q32)*qup3*EL*Mf2(3,Gen4)) + 
     -        3*Alfa*((q12+q22-q32)*qup4*EL*Mf2(3,Gen4)) + 
     -        12*Alfa*(q22*qup2*EL*TB*Mf2(3,Gen4))))/(SB*SW*TB)+
     -        (Alfa*(9 + SW2*(-12 + 16*SW2)))/(12.D0*CW2*MW*Pi*SW2)*
     -        (qup2*EL*Bval(bb0,iintQ3(Gen4))*Mf2(3,Gen4))/(SB*SW)
     -        - (Alfa*(9 + SW2*(-24 + 32*SW2)))/(6.D0*CW2*MW*Pi*SW2)*
     -        (qup2*EL*Cval(cc00,iintQ4(Gen4))*Mf2(3,Gen4))/
     -        (SB*SW) + Cval(cc0,iintQ4(Gen4))*((Alfa*(3 - 4*SW2))/
     -        (12.D0*CW2*MW*Pi)*(((q12+q22-q32)*qup3 + (q12+q22-q32)*
     -        qup4)*EL*Mf2(3,Gen4))/(SB*SW*TB) + (Alfa*(3 - 4*SW2))/
     -        (3.D0*CW2*MW*Pi)*(qup2*EL*Mf2(3,Gen4)**2)/(SB*SW) + 
     -        (Alfa*(9 + SW2*(-12 + 16*SW2)))/(12.D0*CW2*MW*Pi*SW2)*
     -        (qup2*EL*Mf2(3,Gen4)**2)/(SB*SW)) - 
     -        (Alfa*(9 + SW2*(-6 + 4*SW2)))/(12.D0*CW2*MW*Pi*SW2)*
     -        (qup1*EL*Bval(bb0,iintQ1(Gen4))*Mf2(4,Gen4))/(CB*SW)
     -        + (Alfa*(9 + SW2*(-12 + 8*SW2)))/(6.D0*CW2*MW*Pi*SW2)*
     -        (qup1*EL*Cval(cc00,iintQ2(Gen4))*Mf2(4,Gen4))/
     -        (CB*SW) - (3*Alfa)/(32.D0*CW2*MW*Pi*SW2)*
     -        (qup5*EL*Cval(cc2,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW)
         
         loopa1 = loopa1 - (3*Alfa)/(32.D0*CW2*MW*Pi*SW2)*
     -        (qup6*EL*Cval(cc1,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) + 
     -        Cval(cc0,iintQ2(Gen4))*(-((Alfa*(3 - 2*SW2))/(6.D0*CW2*MW*
     -        Pi)*(qup1*EL*Mf2(4,Gen4)**2)/(CB*SW)) - 
     -        (Alfa*(9 + SW2*(-6 + 4*SW2)))/(12.D0*CW2*MW*Pi*SW2)*
     -        (qup1*EL*Mf2(4,Gen4)**2)/(CB*SW) - 
     -        (Alfa*(3 - 2*SW2))/(24.D0*CW2*MW*Pi)*(EL*Mf2(4,Gen4)*
     -        ((q12+q22-q32)*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -        (q12+q22-q32)*(cI*(SB*ZHggs(HiggsType,3)) + 
     -        SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2))))/(CB*SW))


         qup5 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -        ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 2*(cI*(SB*
     -        ZHggs(HiggsType,3)) + SA*
     -        ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))
         qup6 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -        ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 2*(cI*(SB*
     -        ZHggs(HiggsType,3)) + SA*
     -        ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))

         loopa2 = loopa2 + 1/(32.D0*CW2*MW*Pi*SW2)*
     -        (Cval(cc2,iintQ4(Gen4))*(3*Alfa*(2*qup3*EL*Mf2(3,Gen4)) + 
     -        3*Alfa*(2*qup4*EL*Mf2(3,Gen4))))/(SB*SW*TB)+
     -        1/(32.D0*CW2*MW*Pi*SW2)*(Cval(cc1,iintQ4(Gen4))*(3*Alfa*
     -        (2*qup3*EL*Mf2(3,Gen4)) +3*Alfa*(2*qup4*EL*Mf2(3,Gen4))))/
     -        (SB*SW*TB) + (Alfa*(9 + SW2*(-24 + 32*SW2)))/(6.D0*CW2*MW*
     -        Pi*SW2)*(qup2*EL*Cval(cc12,iintQ4(Gen4))*Mf2(3,Gen4))/(SB*
     -        SW) + Cval(cc0,iintQ4(Gen4))*((Alfa*(3 - 4*SW2))/(12.D0*
     -        CW2*MW*Pi)*((2*qup3 +2*qup4)*EL*Mf2(3,Gen4))/(SB*SW*TB)) - 
     -        (Alfa*(9 + SW2*(-12 + 8*SW2)))/(6.D0*CW2*MW*Pi*SW2)*
     -        (qup1*EL*Cval(cc12,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) -
     -        (3*Alfa)/(32.D0*CW2*MW*Pi*SW2)*
     -        (qup5*EL*Cval(cc2,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW)
         
         loopa2 = loopa2 - (3*Alfa)/(32.D0*CW2*MW*Pi*SW2)*
     -        (qup6*EL*Cval(cc1,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) + 
     -        Cval(cc0,iintQ2(Gen4))*(-(Alfa*(3 - 2*SW2))/(24.D0*CW2*MW*
     -        Pi)*(EL*Mf2(4,Gen4)*(2*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -        2*(cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2))))/(CB*SW))


         qup5 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -        ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + (-2)*(cI*
     -        (SB*ZHggs(HiggsType,3)) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))
         qup6 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -        ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + (-2)*(cI*
     -        (SB*ZHggs(HiggsType,3)) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))

         loopa3 = loopa3 + 1/(32.D0*CW2*MW*Pi*SW2)*
     -        (Cval(cc2,iintQ4(Gen4))*(3*Alfa*((-2)*qup3*EL*Mf2(3,Gen4))+ 
     -        3*Alfa*(2*qup4*EL*Mf2(3,Gen4))))/(SB*SW*TB)+
     -        1/(32.D0*CW2*MW*Pi*SW2)*(Cval(cc1,iintQ4(Gen4))*
     -        (3*Alfa*((-2)*qup3*EL*Mf2(3,Gen4)) + 
     -         3*Alfa*(2*qup4*EL*Mf2(3,Gen4))))/(SB*SW*TB) + 
     -         Cval(cc0,iintQ4(Gen4))*
     -        ((Alfa*(3 - 4*SW2))/(12.D0*CW2*MW*Pi)*
     -        (((-2)*qup3 + 2*qup4)*EL*Mf2(3,Gen4))/(SB*SW*TB)) - 
     -        (3*Alfa)/(32.D0*CW2*MW*Pi*SW2)*
     -        (qup5*EL*Cval(cc2,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW)
         
         loopa3 = loopa3 - (3*Alfa)/(32.D0*CW2*MW*Pi*SW2)*
     -        (qup6*EL*Cval(cc1,iintQ2(Gen4))*Mf2(4,Gen4))/(CB*SW) + 
     -        Cval(cc0,iintQ2(Gen4))*
     -        (-(Alfa*(3 - 2*SW2))/(24.D0*CW2*MW*Pi)*(EL*Mf2(4,Gen4)*
     -        (2*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -        (-2)*(cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2))))/(CB*SW))
      end do
      
      do Sfe5 = 1, 2
         do Sfe4 = 1, 2
            do Gen4 = g4ini, 3       
               qup7 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
               qup8 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
               qup9 = -(SAB*ZHggs(HiggsType,1))+CAB*ZHggs(HiggsType,2)
               qup10 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
               qup11 = -(cI*(SB*ZHggs(HiggsType,3)))+SA*
     -              ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
               qup12 = cI*(SB*ZHggs(HiggsType,3)) + SA*
     -              ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)

               loopa1 = loopa1 - Alfa/(144.D0*CW*CW2*MW*Pi*SW2)*
     -              (EL*Bval(bb0,iintQ6(Gen4,Sfe4,Sfe5))*((3- 4*SW2)**2*
     -              (USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)) + 
     -              16*SW2**2*(USf(Sfe4,2,3,Gen4)*USfC(Sfe5,2,3,Gen4)))*
     -              (USf(Sfe5,1,3,Gen4)*(3*CW*(Mf(3,Gen4)*
     -              USfC(Sfe4,2,3,Gen4)*(qup11*MUE*TB +AfC(3,Gen4,Gen4)*
     -              (-(cI*(SB*ZHggs(HiggsType,3))) + CA*TB*
     -              ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2))))+
     -              TB*(MW*MZ*(-3 + 4*SW2)*(qup8*SB) + 
     -              6*CW*(qup10*Mf2(3,Gen4)))*USfC(Sfe4,1,3,Gen4)) +
     -              USf(Sfe5,2,3,Gen4)*(3*CW*(Mf(3,Gen4)*
     -              USfC(Sfe4,1,3,Gen4)*(qup12*MUEC*TB +Af(3,Gen4,Gen4)*
     -              (cI*(SB*ZHggs(HiggsType,3)) + CA*TB*
     -              ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)))) + 
     -              2*TB*(2*MW*MZ*SW2*(qup9*SB) + 3*CW*(qup10*
     -              Mf2(3,Gen4)))*USfC(Sfe4,2,3,Gen4))))/(SB*SW*TB) + 
     -              Alfa/(144.D0*CW*CW2*MW*Pi*SW2)*(EL*
     -              Bval(bb0,iintQ5(Gen4,Sfe4,Sfe5))*((3 - 2*SW2)**2*
     -              (USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)) + 
     -              4*SW2**2*(USf(Sfe4,2,4,Gen4)*USfC(Sfe5,2,4,Gen4)))*
     -              (USf(Sfe5,1,4,Gen4)*(3*CW*(Mf(4,Gen4)*
     -              USfC(Sfe4,2,4,Gen4)*(qup12*AfC(4,Gen4,Gen4) + MUE*
     -              (cI*(CB*ZHggs(HiggsType,3))+CA*ZHggs(HiggsType,1)+ 
     -              SA*ZHggs(HiggsType,2)))) + (MW*MZ*(-3 + 2*SW2)*(CB*
     -              qup8)+6*CW*(qup7*Mf2(4,Gen4)))*USfC(Sfe4,1,4,Gen4))+
     -              USf(Sfe5,2,4,Gen4)*(3*CW*(Mf(4,Gen4)*
     -              USfC(Sfe4,1,4,Gen4)*(qup11*Af(4,Gen4,Gen4) + MUEC*
     -              (-(cI*(CB*ZHggs(HiggsType,3))) + CA*
     -              ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)))) + 2*
     -              (MW*MZ*SW2*(CB*qup9) + 3*CW*(qup7*Mf2(4,Gen4)))*
     -              USfC(Sfe4,2,4,Gen4))))/(CB*SW)
            end do 
         end do
      end do

      do Sfe6 = 1,2
         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = g4ini,3
                  qup13 = SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                  qup14 = SAB*ZHggs(HiggsType,1) - CAB*
     -                 ZHggs(HiggsType,2)
                  qup15 = -(SAB*ZHggs(HiggsType,1)) + CAB*
     -                 ZHggs(HiggsType,2)
                  qup16 = CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                  qup17 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -                 CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                  qup18 = cI*(CB*ZHggs(HiggsType,3))+CA*
     -                 ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                  qup19 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -                 SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                  qup20 = -(cI*(SB*ZHggs(HiggsType,3))) + CA*TB*
     -                 ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)
                  qup21 = cI*(SB*ZHggs(HiggsType,3))+SA*
     -                 ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
                  qup22 = cI*(SB*ZHggs(HiggsType,3)) + CA*TB*
     -                 ZHggs(HiggsType,1) + SA*TB*ZHggs(HiggsType,2)
                  qup23 = MW*MZ*(-3 + 4*SW2)*(qup14*SB) + 
     -                 6*CW*(qup16*Mf2(3,Gen4))
                  qup24 = MW*MZ*(-3 + 2*SW2)*(CB*qup14) + 
     -                 6*CW*(qup13*Mf2(4,Gen4))
                  qup25 = SW2*(4*USf(Sfe5,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)+ 
     -                 4*USf(Sfe5,2,3,Gen4)*USfC(Sfe4,2,3,Gen4)) - 
     -                 3*USf(Sfe5,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)
                  qup26 = SW2*(4*USf(Sfe6,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)+ 
     -                 4*USf(Sfe6,2,3,Gen4)*USfC(Sfe4,2,3,Gen4)) - 
     -                 3*USf(Sfe6,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)
                  qup27 = SW2*(2*USf(Sfe5,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)+ 
     -                 2*USf(Sfe5,2,4,Gen4)*USfC(Sfe4,2,4,Gen4)) - 
     -                 3*USf(Sfe5,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)
                  qup28 = SW2*(2*USf(Sfe6,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)+ 
     -                 2*USf(Sfe6,2,4,Gen4)*USfC(Sfe4,2,4,Gen4)) - 
     -                 3*USf(Sfe6,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)
                  qup29 = SW2*(4*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)+ 
     -                 4*USf(Sfe4,2,3,Gen4)*USfC(Sfe5,2,3,Gen4)) - 
     -                 3*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)
                  qup30 = SW2*(2*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)+ 
     -                 2*USf(Sfe4,2,4,Gen4)*USfC(Sfe5,2,4,Gen4)) - 
     -                 3*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)
                  qup31 = SW2*(4*USf(Sfe4,1,3,Gen4)*USfC(Sfe6,1,3,Gen4)+ 
     -                 4*USf(Sfe4,2,3,Gen4)*USfC(Sfe6,2,3,Gen4)) - 
     -                 3*USf(Sfe4,1,3,Gen4)*USfC(Sfe6,1,3,Gen4)
                  qup32 = SW2*(2*USf(Sfe4,1,4,Gen4)*USfC(Sfe6,1,4,Gen4)+ 
     -                 2*USf(Sfe4,2,4,Gen4)*USfC(Sfe6,2,4,Gen4)) - 
     -                 3*USf(Sfe4,1,4,Gen4)*USfC(Sfe6,1,4,Gen4)
                  qup33 = 3*CW*((qup19*MUE*TB + qup20*AfC(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe5,2,3,Gen4)) + 
     -                 qup23*TB*USfC(Sfe5,1,3,Gen4)
                  qup34 = 3*CW*((qup21*MUEC*TB + qup22*Af(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe5,1,3,Gen4)) + 2*TB*(2*MW*MZ*
     -                 SW2*(qup15*SB) + 3*CW*(qup16*Mf2(3,Gen4)))*
     -                 USfC(Sfe5,2,3,Gen4)
                  qup35 = 3*CW*((qup18*MUE + qup21*AfC(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe5,2,4,Gen4)) + 
     -                 qup24*USfC(Sfe5,1,4,Gen4)
                  qup36 = 3*CW*((qup17*MUEC + qup19*Af(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe5,1,4,Gen4)) + 2*(MW*MZ*SW2*
     -                 (CB*qup15) + 3*CW*(qup13*Mf2(4,Gen4)))*
     -                 USfC(Sfe5,2,4,Gen4)
                  qup37 = 3*CW*((qup19*MUE*TB + qup20*AfC(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe6,2,3,Gen4)) + 
     -                 qup23*TB*USfC(Sfe6,1,3,Gen4)
                  qup38 = 3*CW*((qup21*MUEC*TB + qup22*Af(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe6,1,3,Gen4)) + 2*TB*(2*MW*MZ*
     -                 SW2*(qup15*SB) + 3*CW*(qup16*Mf2(3,Gen4)))*
     -                 USfC(Sfe6,2,3,Gen4)
                  qup39 = 3*CW*((qup18*MUE + qup21*AfC(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe6,2,4,Gen4)) + 
     -                 qup24*USfC(Sfe6,1,4,Gen4)
                  qup40 = 3*CW*((qup17*MUEC + qup19*Af(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe6,1,4,Gen4)) + 2*(MW*MZ*SW2*
     -                 (CB*qup15) + 3*CW*(qup13*Mf2(4,Gen4)))*
     -                 USfC(Sfe6,2,4,Gen4)

                  loopa1 = loopa1 + Alfa/(72.D0*CW*CW2*MW*Pi*SW2)*
     -                 (EL*Cval(cc00,iintQ8(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 (qup26*qup29*(qup37*USf(Sfe5,1,3,Gen4) + 
     -                 qup38*USf(Sfe5,2,3,Gen4)) + qup25*qup31*
     -                 (qup33*USf(Sfe6,1,3,Gen4) + 
     -                 qup34*USf(Sfe6,2,3,Gen4))))/(SB*SW*TB) - 
     -                 Alfa/(72.D0*CW*CW2*MW*Pi*SW2)*
     -                 (EL*Cval(cc00,iintQ7(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 (qup28*qup30*(qup39*USf(Sfe5,1,4,Gen4) + 
     -                 qup40*USf(Sfe5,2,4,Gen4)) + qup27*qup32*
     -                 (qup35*USf(Sfe6,1,4,Gen4) + 
     -                 qup36*USf(Sfe6,2,4,Gen4))))/(CB*SW) 


                  loopa2 = loopa2 - Alfa/(72.D0*CW*CW2*MW*Pi*SW2)*
     -                 (EL*Cval(cc12,iintQ8(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 (qup26*qup29*(qup37*USf(Sfe5,1,3,Gen4) + 
     -                 qup38*USf(Sfe5,2,3,Gen4)) + qup25*qup31*
     -                 (qup33*USf(Sfe6,1,3,Gen4) + 
     -                 qup34*USf(Sfe6,2,3,Gen4))))/(SB*SW*TB) + 
     -                 Alfa/(72.D0*CW*CW2*MW*Pi*SW2)*
     -                 (EL*Cval(cc12,iintQ7(Gen4,Sfe4,Sfe5,Sfe6))*
     -                 (qup28*qup30*(qup39*USf(Sfe5,1,4,Gen4) + 
     -                 qup40*USf(Sfe5,2,4,Gen4)) + qup27*qup32*
     -                 (qup35*USf(Sfe6,1,4,Gen4) + 
     -                 qup36*USf(Sfe6,2,4,Gen4))))/(CB*SW)
               end do
            end do
         end do
      end do


** The lepton and slepton contributions
      if (sector .ne. 1) then
         do Gen4 = 1, 3
            lup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
            lup2 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)

            lup3 = q12*(4*SA*ZHggs(HiggsType,1) - 4*CA*
     -           ZHggs(HiggsType,2)) + 
     -           (q12+q22-q32)*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -           (q12+q22-q32)*(cI*(SB*ZHggs(HiggsType,3)) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))
            lup4 = q22*(4*SA*ZHggs(HiggsType,1) - 4*CA*
     -           ZHggs(HiggsType,2)) + 
     -           (q12+q22-q32)*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -           (q12+q22-q32)*(cI*(SB*ZHggs(HiggsType,3)) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))

            loopa1 = loopa1 + (Alfa*MZ)/(16.D0*CW*CW2*Pi*SW2)*
     -           (lup2*EL*Bval(bb0,iintL3(Gen4)))/SW - 
     -           (Alfa*MZ)/(4.D0*CW*CW2*Pi*SW2)*
     -           (lup2*EL*Cval(cc00,iintL4(Gen4)))/SW - 
     -           (Alfa*(1 + SW2*(-2 + 4*SW2)))/(4.D0*CW2*MW*Pi*SW2)*
     -           (lup1*EL*Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4))/(CB*SW) +
     -           (Alfa*(1 + SW2*(-4 + 8*SW2)))/(2.D0*CW2*MW*Pi*SW2)*
     -           (lup1*EL*Cval(cc00,iintL2(Gen4))*Mf2(2,Gen4))/
     -           (CB*SW) - Alfa/(32.D0*CW2*MW*Pi*SW2)*
     -           (lup4*EL*Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) - 
     -           Alfa/(32.D0*CW2*MW*Pi*SW2)*
     -           (lup3*EL*Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) + 
     -           Cval(cc0,iintL2(Gen4))*
     -           (-((Alfa*(1 - 2*SW2))/(2.D0*CW2*MW*Pi)*
     -           (lup1*EL*Mf2(2,Gen4)**2)/(CB*SW)) - 
     -           (Alfa*(1 + SW2*(-2 + 4*SW2)))/(4.D0*CW2*MW*Pi*SW2)*
     -           (lup1*EL*Mf2(2,Gen4)**2)/(CB*SW) - 
     -           (Alfa*(1 - 2*SW2))/(8.D0*CW2*MW*Pi)*(EL*Mf2(2,Gen4)*
     -           ((q12+q22-q32)*(-(cI*(SB*ZHggs(HiggsType,3))) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -           (q12+q22-q32)*(cI*(SB*ZHggs(HiggsType,3)) + SA*
     -           ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2))))/(CB*SW))


            lup3 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -           ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 2*(cI*
     -           (SB*ZHggs(HiggsType,3)) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))            
            lup4 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) +SA*
     -           ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 2*(cI*
     -           (SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - CA*
     -           ZHggs(HiggsType,2))

            loopa2 = loopa2 + (Alfa*MZ)/(4.D0*CW*CW2*Pi*SW2)*(lup2*EL*
     -           Cval(cc12,iintL4(Gen4)))/SW-Alfa/(32.D0*CW2*MW*Pi*SW2)*
     -           (lup4*EL*Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) - 
     -           (Alfa*(1 + SW2*(-4 + 8*SW2)))/(2.D0*CW2*MW*Pi*SW2)*
     -           (lup1*EL*Cval(cc12,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) -
     -           Alfa/(32.D0*CW2*MW*Pi*SW2)*(lup3*EL*
     -           Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) + 
     -           Cval(cc0,iintL2(Gen4))*(-(Alfa*(1 - 2*SW2))/(8.D0*CW2*
     -           MW*Pi)*(EL*Mf2(2,Gen4)*(2*(-(cI*(SB*
     -           ZHggs(HiggsType,3)))+ SA*ZHggs(HiggsType,1) - CA*
     -           ZHggs(HiggsType,2)) + 2*(cI*(SB*ZHggs(HiggsType,3)) +
     -           SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2))))/(CB*SW))


            lup3 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) +SA*
     -           ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + (-2)*
     -           (cI*(SB*ZHggs(HiggsType,3)) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))
            lup4 = 2*(-(cI*(SB*ZHggs(HiggsType,3))) +SA*
     -           ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + (-2)*
     -           (cI*(SB*ZHggs(HiggsType,3)) + 
     -           SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2))

            loopa3 = loopa3 - Alfa/(32.D0*CW2*MW*Pi*SW2)*
     -           (lup4*EL*Cval(cc1,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) - 
     -           Alfa/(32.D0*CW2*MW*Pi*SW2)*(lup3*EL*
     -           Cval(cc2,iintL2(Gen4))*Mf2(2,Gen4))/(CB*SW) + 
     -           Cval(cc0,iintL2(Gen4))*(-(Alfa*(1 - 2*SW2))/(8.D0*CW2*
     -           MW*Pi)*(EL*Mf2(2,Gen4)*(2*(-(cI*(SB*
     -           ZHggs(HiggsType,3)))+ SA*ZHggs(HiggsType,1) - CA*
     -           ZHggs(HiggsType,2)) + (-2)*(cI*(SB*
     -           ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2))))/(CB*SW))
         enddo
            
         do Sfe5 = 1, 2
            do Sfe4 = 1, 2
               do Gen4 = 1, 3
                  lup5 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)

                  loopa1 = loopa1 + Alfa/(16.D0*CW*CW2*MW*Pi*SW2)*
     -                 (EL*Bval(bb0,iintL5(Gen4,Sfe4,Sfe5))*
     -                 ((1 - 2*SW2)**2*(USf(Sfe4,1,2,Gen4)*
     -                 USfC(Sfe5,1,2,Gen4)) + 4*SW2**2*
     -                 (USf(Sfe4,2,2,Gen4)*USfC(Sfe5,2,2,Gen4)))*
     -                 (USf(Sfe5,1,2,Gen4)*(CW*(Mf(2,Gen4)*
     -                 USfC(Sfe4,2,2,Gen4)*(AfC(2,Gen4,Gen4)*(cI*(SB*
     -                 ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -                 CA*ZHggs(HiggsType,2)) + MUE*(cI*(CB*
     -                 ZHggs(HiggsType,3)) + CA*ZHggs(HiggsType,1) + 
     -                 SA*ZHggs(HiggsType,2))))+(2*CW*(lup5*
     -                 Mf2(2,Gen4)) + MW*MZ*(-1 + 2*SW2)*(CB*(SAB*
     -                 ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2))))*
     -                 USfC(Sfe4,1,2,Gen4)) + USf(Sfe5,2,2,Gen4)*(CW*
     -                 (Mf(2,Gen4)*USfC(Sfe4,1,2,Gen4)*(Af(2,Gen4,Gen4)*
     -                 (-(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -                 ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -                 MUEC*(-(cI*(CB*ZHggs(HiggsType,3))) + CA*
     -                 ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))))+
     -                 2*(CW*(lup5*Mf2(2,Gen4)) + MW*MZ*SW2*(CB*(-(SAB*
     -                 ZHggs(HiggsType,1)) + CAB*ZHggs(HiggsType,2))))*
     -                 USfC(Sfe4,2,2,Gen4))))/(CB*SW)
               enddo
            enddo
         enddo

            
         do Sfe6 = 1, 2
            do Sfe5 = 1, 2
               do Sfe4 = 1, 2
                  do Gen4 = 1, 3
                     lup6 = SA*ZHggs(HiggsType,1) - CA*
     -                    ZHggs(HiggsType,2)
                     lup7 = SAB*ZHggs(HiggsType,1) - CAB*
     -                    ZHggs(HiggsType,2)
                     lup8 = -(SAB*ZHggs(HiggsType,1)) + 
     -                    CAB*ZHggs(HiggsType,2)
                     lup9 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -                    CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                     lup10 = cI*(CB*ZHggs(HiggsType,3)) + 
     -                    CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                     lup11 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -                    SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                     lup12 = cI*(SB*ZHggs(HiggsType,3)) + 
     -                    SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                     lup13 = MW*MZ*(-1 + 2*SW2)*(CB*lup7) + 
     -                    2*CW*(lup6*Mf2(2,Gen4))
                     lup14 = SW2*(2*USf(Sfe5,1,2,Gen4)*
     -                    USfC(Sfe4,1,2,Gen4) + 
     -                    2*USf(Sfe5,2,2,Gen4)*USfC(Sfe4,2,2,Gen4)) - 
     -                    USf(Sfe5,1,2,Gen4)*USfC(Sfe4,1,2,Gen4)
                     lup15 = SW2*(2*USf(Sfe6,1,2,Gen4)*
     -                    USfC(Sfe4,1,2,Gen4) + 
     -                    2*USf(Sfe6,2,2,Gen4)*USfC(Sfe4,2,2,Gen4)) - 
     -                    USf(Sfe6,1,2,Gen4)*USfC(Sfe4,1,2,Gen4)
                     lup16 = SW2*(2*USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4) + 
     -                    2*USf(Sfe4,2,2,Gen4)*USfC(Sfe5,2,2,Gen4)) - 
     -                    USf(Sfe4,1,2,Gen4)*USfC(Sfe5,1,2,Gen4)
                     lup17 = SW2*(2*USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe6,1,2,Gen4) + 
     -                    2*USf(Sfe4,2,2,Gen4)*USfC(Sfe6,2,2,Gen4)) - 
     -                    USf(Sfe4,1,2,Gen4)*USfC(Sfe6,1,2,Gen4)
                     lup18 = CW*((lup10*MUE + lup12*AfC(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe5,2,2,Gen4)) + 
     -                    lup13*USfC(Sfe5,1,2,Gen4)
                     lup19 = CW*((lup9*MUEC + lup11*Af(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe5,1,2,Gen4)) + 2*(MW*MZ*
     -                    SW2*(CB*lup8) + CW*(lup6*Mf2(2,Gen4)))*
     -                    USfC(Sfe5,2,2,Gen4)
                     lup20 = CW*((lup10*MUE + lup12*AfC(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe6,2,2,Gen4)) + 
     -                    lup13*USfC(Sfe6,1,2,Gen4)
                     lup21 = CW*((lup9*MUEC + lup11*Af(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe6,1,2,Gen4)) + 
     -                    2*(MW*MZ*SW2*(CB*lup8) + 
     -                    CW*(lup6*Mf2(2,Gen4)))*USfC(Sfe6,2,2,Gen4)
                     
                     loopa1 = loopa1 - Alfa/(8.D0*CW*CW2*MW*Pi*SW2)*
     -                    (EL*Cval(cc00,iintL6(Gen4,Sfe4,Sfe5,Sfe6))*
     -                    (lup15*lup16*(lup20*USf(Sfe5,1,2,Gen4) + 
     -                    lup21*USf(Sfe5,2,2,Gen4)) + lup14*lup17*
     -                    (lup18*USf(Sfe6,1,2,Gen4) + 
     -                    lup19*USf(Sfe6,2,2,Gen4))))/(CB*SW) 


                     loopa2 = loopa2 + Alfa/(8.D0*CW*CW2*MW*Pi*SW2)*
     -                    (EL*Cval(cc12,iintL6(Gen4,Sfe4,Sfe5,Sfe6))*
     -                    (lup15*lup16*(lup20*USf(Sfe5,1,2,Gen4) + 
     -                    lup21*USf(Sfe5,2,2,Gen4)) + lup14*lup17*
     -                    (lup18*USf(Sfe6,1,2,Gen4) + 
     -                    lup19*USf(Sfe6,2,2,Gen4))))/(CB*SW)
                  enddo
               enddo
            enddo
         enddo

      end if


** fermion values of the formfactors
      a2vtxZfm = -loopa2
      a1vtxZfm = loopa1 + (loopa2*q1q2) + ctFM
      a3vtxZfm = loopa3*(0.D0,-1.D0)


** And the SM bosons and their superpartners ...
      if (sector .ge. 3) then
         bup1 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
         bup2 = SBA*ZHggs(HiggsType,1) + CBA*ZHggs(HiggsType,2)
         bup3 = 2*CW2*(SBA*ZHggs(HiggsType,1)) + 
     -        2*CW2*(CBA*ZHggs(HiggsType,2)) + 
     -        C2B*(SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2))

         loopa1 = loopa1 + Bval(bb0,iintB6)*
     -        ((Alfa*MW*(CW2 - SW2)**2)/(16.D0*CW2**2*Pi*SW2)*
     -        (C2B*bup1*EL)/SW + (Alfa*MW)/(2.D0*Pi)*(bup2*EL)/SW - 
     -        (Alfa*CW2*MW)/(2.D0*Pi*SW2)*(bup2*EL)/SW) + 
     -        Cval(cc00,iintB7)*(-((Alfa*MW*(CW2 - SW2)**2)/
     -        (4.D0*CW2**2*Pi*SW2)*(C2B*bup1*EL)/SW) - 
     -        (Alfa*MW)/(2.D0*Pi)*(bup2*EL)/SW + 
     -        (9*Alfa*CW2*MW)/(2.D0*Pi*SW2)*(bup2*EL)/SW + 
     -        (Alfa*MW*(-CW2 + SW2))/(2.D0*CW2*Pi)*(bup2*EL)/SW)+
     -        Cval(cc0,iintB7)*(-((Alfa*MW)/(4.D0*Pi)*
     -        (((q12-q32) + (q22-q32))*bup2*EL)/SW) + 
     -        (5*Alfa*CW2*MW)/(4.D0*Pi*SW2)*((q12+q22-q32)*bup2*EL)/SW + 
     -        (Alfa*MW*MW2*SW2)/(4.D0*CW2**2*Pi)*(C2B*bup1*EL)/SW + 
     -        (Alfa*MW*MW2)/(2.D0*Pi)*(bup2*EL)/SW - 
     -        (Alfa*MW*MW2*(-2*CW2**2 + SW2**2))/(2.D0*CW2*Pi*SW2)*
     -        (bup2*EL)/SW) + (Alfa*MW)/(32.D0*CW2**2*Pi*SW2)*
     -        (C2B*bup1*EL*Bval(bb0,iintB1))/SW - 
     -        (Alfa*MW)/(4.D0*CW2*Pi)*
     -        (bup2*EL*SW*(Bval(bb0,iintB2) + Bval(bb0,iintB3))) - 
     -        (Alfa*MW*(CW2 - SW2)**2)/(16.D0*CW2**2*Pi*SW2)*
     -        (bup3*EL*Bval(bb0,iintB4))/SW + 
     -        (Alfa*MW*(CW2 - SW2)**2)/(4.D0*CW2**2*Pi*SW2)*
     -        (bup3*EL*Cval(cc00,iintB5))/SW + 1/(4.D0*CW2*Pi*SW2)*
     -        (Cval(cc1,iintB7)*(Alfa*CW2*MW*SW2*((q12-q22-q32)*bup2*EL)- 
     -        Alfa*CW2*MW*SW2*((q12-q22-q32)*bup2*EL) + 
     -        Alfa*CW2**2*MW*((q12+3*q22-q32)*bup2*EL)))/SW + 
     -        1/(4.D0*CW2*Pi*SW2)*(Cval(cc2,iintB7)*
     -        (Alfa*CW2*MW*SW2*((q12-q22+q32)*bup2*EL) - 
     -        Alfa*CW2*MW*SW2*((q12-q22+q32)*bup2*EL) + 3*Alfa*CW2**2*
     -        MW*((q12+1/3.D0*q22-1/3.D0*q32)*bup2*EL)))/SW


         loopa2 = loopa2 + Cval(cc12,iintB7)*((Alfa*MW*(CW2 - SW2)**2)/
     -        (4.D0*CW2**2*Pi*SW2)*(C2B*bup1*EL)/SW + 
     -        (Alfa*MW)/(2.D0*Pi)*(bup2*EL)/SW - 
     -        (9*Alfa*CW2*MW)/(2.D0*Pi*SW2)*(bup2*EL)/SW - 
     -        (Alfa*MW*(-CW2 + SW2))/(2.D0*CW2*Pi)*(bup2*EL)/SW) +
     -        Cval(cc0,iintB7)*(-((Alfa*MW)/(4.D0*Pi)*(8*bup2*EL)/SW) + 
     -        (5*Alfa*CW2*MW)/(4.D0*Pi*SW2)*((8/5D0)*bup2*EL)/SW) - 
     -        (Alfa*MW*(CW2 - SW2)**2)/(4.D0*CW2**2*Pi*SW2)*
     -        (bup3*EL*Cval(cc12,iintB5))/SW + 1/(4.D0*CW2*Pi*SW2)*
     -        (Cval(cc1,iintB7)*(Alfa*CW2*MW*SW2*((-2)*bup2*EL) - 
     -        Alfa*CW2*MW*SW2*(4*bup2*EL) + Alfa*CW2**2*MW*((-2)*bup2*
     -        EL) + Alfa*((2*CW2*MW*SW2*(bup2*EL) - 2*MW*SW2**2*(bup2*
     -        EL)))))/SW + 1/(4.D0*CW2*Pi*SW2)*(Cval(cc2,iintB7)*
     -        (Alfa*CW2*MW*SW2*((-4)*bup2*EL) - Alfa*CW2*MW*SW2*(2*bup2*
     -        EL) + 3*Alfa*CW2**2*MW*((-2/3D0)*bup2*EL) + Alfa*((2*CW2*
     -        MW*SW2*(bup2*EL) - 2*MW*SW2**2*(bup2*EL)))))/SW

         
** Higgs contributions ..
         do Hig4 = 1,3
            bup18 = SBA*UHggs(Hig4,1) + CBA*UHggs(Hig4,2)
            bup19 = SBA*ZHggs(HiggsType,1) + CBA*ZHggs(HiggsType,2)
            bup20 = 2*SBA*ZHggs(HiggsType,1) + 2*CBA*
     -           ZHggs(HiggsType,2) + C2B*(-(SAB*ZHggs(HiggsType,1)) + 
     -           CAB*ZHggs(HiggsType,2))

            loopa1 = loopa1 + (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -           (bup18**2*bup20*EL*Cval(cc00,iintB10(Hig4)))/SW + 
     -           (Alfa*MW)/(4.D0*CW2**3*Pi*SW2)*(bup18**2*bup19*EL*
     -           Cval(cc0,iintB10(Hig4))*(-MW2))/SW - (Alfa*MW)/(8.D0*
     -           CW2**2*Pi*SW2)*(bup18*EL*(Bval(bb0,iintB11(Hig4)) + 
     -           Bval(bb0,iintB12(Hig4)))*(UHggs(Hig4,1)*
     -           ZHggs(HiggsType,1)+UHggs(Hig4,2)*ZHggs(HiggsType,2)+ 
     -           UHggs(Hig4,3)*ZHggs(HiggsType,3)))/SW


            loopa2 = loopa2 + (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -           (bup18**2*bup20*EL*(Cval(cc1,iintB10(Hig4)) + 
     -           Cval(cc11,iintB10(Hig4)) + Cval(cc12,iintB10(Hig4))))/
     -           SW + (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -           (bup18**2*bup19*EL*Cval(cc2,iintB10(Hig4)))/SW + 
     -           (Alfa*MW)/(4.D0*CW2**3*Pi*SW2)*(bup18**2*bup19*EL*
     -           Cval(cc0,iintB10(Hig4))*CW2)/SW 
         enddo


         do Hig5 = 1,3
            do Hig4 = 1,3
               bup21 = SBA*UHggs(Hig4,1) + CBA*UHggs(Hig4,2)
               bup22 = CBA*UHggs(Hig4,1) - SBA*UHggs(Hig4,2)
               bup23 = SAB*UHggs(Hig5,1) - CAB*UHggs(Hig5,2)
               bup24 = SBA*UHggs(Hig5,1) + CBA*UHggs(Hig5,2)
               bup25 = -(CBA*UHggs(Hig5,1)) + SBA*UHggs(Hig5,2)
               bup26 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
               bup27 = 2*CAB*S2A*ZHggs(HiggsType,1) + C2A*(SAB*
     -              ZHggs(HiggsType,1) - 3*CAB*ZHggs(HiggsType,2))
               bup28 = -2*S2A*SAB*ZHggs(HiggsType,1) + 2*CAB*S2A*
     -              ZHggs(HiggsType,2) + C2A*(CAB*ZHggs(HiggsType,1) + 
     -              SAB*ZHggs(HiggsType,2))
               bup29 = 2*S2A*SAB*ZHggs(HiggsType,2) + C2A*(3*SAB*
     -              ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2))
               bup30 = -(CBA*ZHggs(HiggsType,1)) - S2B*SAB*
     -              ZHggs(HiggsType,1) + 
     -              CAB*S2B*ZHggs(HiggsType,2) +SBA*ZHggs(HiggsType,2)
               bup31 = CBA*UHggs(Hig5,1) - S2B*SAB*UHggs(Hig5,1) + 
     -              CAB*S2B*UHggs(Hig5,2) - SBA*UHggs(Hig5,2)
               bup32 = C2B*UHggs(Hig4,3)*(bup26*UHggs(Hig5,3) + bup23*
     -              ZHggs(HiggsType,3)) - UHggs(Hig4,2)*(bup28*
     -              UHggs(Hig5,1) + bup27*UHggs(Hig5,2) + C2B*CAB*
     -              UHggs(Hig5,3)*ZHggs(HiggsType,3))+ UHggs(Hig4,1)*
     -              (bup29*UHggs(Hig5,1) - bup28*UHggs(Hig5,2) + 
     -              C2B*SAB*UHggs(Hig5,3)*ZHggs(HiggsType,3))

               loopa1 = loopa1 - (Alfa*MW*MW2)/(8.D0*CW2**3*Pi*SW2)*
     -              (bup21*bup24*bup32*EL*Cval(cc0,iintB15(Hig4,Hig5)))/
     -              SW + (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(bup21*bup24*
     -              bup32*EL*Cval(cc00,iintB15(Hig4,Hig5)))/SW - (Alfa*
     -              MW)/(32.D0*CW2**2*Pi*SW2)*
     -              (bup32*EL*Bval(bb0,iintB16(Hig4,Hig5))*
     -              (UHggs(Hig4,1)*UHggs(Hig5,1) + 
     -              UHggs(Hig4,2)*UHggs(Hig5,2) + 
     -              UHggs(Hig4,3)*UHggs(Hig5,3)))/SW - 
     -              (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -              (bup21*EL*(Cval(cc00,iintB13(Hig4,Hig5)) + 
     -              Cval(cc00,iintB14(Hig4,Hig5)))*
     -              (bup25*UHggs(Hig4,3) + bup22*UHggs(Hig5,3))*
     -              (bup30*UHggs(Hig5,3)+bup31*ZHggs(HiggsType,3)))/SW 


               loopa2 = loopa2 - (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(bup21*
     -              bup24*bup32*EL*Cval(cc12,iintB15(Hig4,Hig5)))/
     -              SW + (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -              (bup21*EL*(Cval(cc0,iintB13(Hig4,Hig5)) + 
     -              Cval(cc0,iintB14(Hig4,Hig5)) + 
     -              Cval(cc1,iintB13(Hig4,Hig5)) + 
     -              Cval(cc1,iintB14(Hig4,Hig5)))*
     -              (bup25*UHggs(Hig4,3) + bup22*UHggs(Hig5,3))*
     -              (UHggs(Hig5,3)*(CBA*ZHggs(HiggsType,1) - SBA*
     -              ZHggs(HiggsType,2)) + bup25*ZHggs(HiggsType,3)))/
     -              SW - (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -              (bup21*EL*(Cval(cc12,iintB13(Hig4,Hig5)) + 
     -              Cval(cc12,iintB14(Hig4,Hig5)) + 
     -              Cval(cc22,iintB13(Hig4,Hig5)) + 
     -              Cval(cc22,iintB14(Hig4,Hig5)))*
     -              (bup25*UHggs(Hig4,3) +bup22*UHggs(Hig5,3))*(bup30*
     -              UHggs(Hig5,3) + bup31*ZHggs(HiggsType,3)))/SW - 
     -              (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -              (bup21*EL*(Cval(cc2,iintB13(Hig4,Hig5)) + 
     -              Cval(cc2,iintB14(Hig4,Hig5)))*(bup25*UHggs(Hig4,3)+ 
     -              bup22*UHggs(Hig5,3))*(UHggs(Hig5,3)*(-3*CBA*
     -              ZHggs(HiggsType,1) - S2B*SAB*ZHggs(HiggsType,1) + 
     -              CAB*S2B*ZHggs(HiggsType,2) + 3*SBA*
     -              ZHggs(HiggsType,2)) + (3*CBA*UHggs(Hig5,1) - S2B*
     -              SAB*UHggs(Hig5,1) + CAB*S2B*UHggs(Hig5,2) - 3*SBA*
     -              UHggs(Hig5,2))*ZHggs(HiggsType,3)))/SW
            enddo
         enddo

         do Hig6 = 1,3
            do Hig5 = 1,3
               do Hig4 = 1,3
                  bup33 = CBA*UHggs(Hig4,1) - SBA*UHggs(Hig4,2)
                  bup34 = -(CBA*UHggs(Hig5,1)) + SBA*UHggs(Hig5,2)
                  bup35 = SAB*UHggs(Hig6,1) - CAB*UHggs(Hig6,2)
                  bup36 = -(CBA*UHggs(Hig6,1)) + SBA*UHggs(Hig6,2)
                  bup37 = SAB*ZHggs(HiggsType,1) - CAB*
     -                 ZHggs(HiggsType,2)
                  bup38 = 2*CAB*S2A*ZHggs(HiggsType,1) + C2A*(SAB*
     -                 ZHggs(HiggsType,1)-3*CAB*ZHggs(HiggsType,2))
                  bup39 = -2*S2A*SAB*ZHggs(HiggsType,1) + 2*CAB*S2A*
     -                 ZHggs(HiggsType,2) + C2A*(CAB*
     -                 ZHggs(HiggsType,1) + SAB*ZHggs(HiggsType,2))
                  bup40 = 2*S2A*SAB*ZHggs(HiggsType,2) + C2A*(3*SAB*
     -                 ZHggs(HiggsType,1)-CAB*ZHggs(HiggsType,2))
                  bup41 = C2B*UHggs(Hig5,3)*(bup37*UHggs(Hig6,3) + 
     -                 bup35*ZHggs(HiggsType,3)) -UHggs(Hig5,2)*(bup39*
     -                 UHggs(Hig6,1) + bup38*UHggs(Hig6,2) + 
     -                 C2B*CAB*UHggs(Hig6,3)*ZHggs(HiggsType,3)) + 
     -                 UHggs(Hig5,1)*(bup40*UHggs(Hig6,1) - 
     -                 bup39*UHggs(Hig6,2) + 
     -                 C2B*SAB*UHggs(Hig6,3)*ZHggs(HiggsType,3))

                  loopa1 = loopa1 + (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -                 (bup41*EL*Cval(cc00,iintB17(Hig4,Hig5,Hig6))*
     -                 (bup34*UHggs(Hig4,3) + bup33*UHggs(Hig5,3))*
     -                 (bup36*UHggs(Hig4,3) + bup33*UHggs(Hig6,3)))/SW 


                  loopa2 = loopa2 - (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -                 (bup41*EL*Cval(cc12,iintB17(Hig4,Hig5,Hig6))*
     -                 (bup34*UHggs(Hig4,3) + bup33*UHggs(Hig5,3))*
     -                 (bup36*UHggs(Hig4,3) + bup33*UHggs(Hig6,3)))/SW
               enddo
            enddo
         enddo


** chargino/neutralino contributions
         if (sector .ge. 4) then
            do Cha6 = 1,2
               do Cha5 = 1,2
                  do Cha4 = 1,2
                     bup42 = 2*SW2*IndexDelta(Cha4,Cha5)-2*UCha(Cha5,1)*
     -                    UChaC(Cha4,1) - UCha(Cha5,2)*UChaC(Cha4,2)
                     bup43 = 2*SW2*IndexDelta(Cha4,Cha6) - 
     -                    2*UCha(Cha6,1)*UChaC(Cha4,1) - UCha(Cha6,2)*
     -                    UChaC(Cha4,2)
                     bup44 = 2*SW2*IndexDelta(Cha4,Cha5) - 
     -                    2*UCha(Cha4,1)*UChaC(Cha5,1) - UCha(Cha4,2)*
     -                    UChaC(Cha5,2)
                     bup45 = 2*SW2*IndexDelta(Cha4,Cha6) - 
     -                    2*UCha(Cha4,1)*UChaC(Cha6,1) - UCha(Cha4,2)*
     -                    UChaC(Cha6,2)
                     bup46 = 2*SW2*IndexDelta(Cha4,Cha5) - 
     -                    2*VCha(Cha5,1)*VChaC(Cha4,1) - VCha(Cha5,2)*
     -                    VChaC(Cha4,2)
                     bup47 = 2*SW2*IndexDelta(Cha4,Cha6) - 
     -                    2*VCha(Cha6,1)*VChaC(Cha4,1) - VCha(Cha6,2)*
     -                    VChaC(Cha4,2)
                     bup48 = 2*SW2*IndexDelta(Cha4,Cha5) - 
     -                    2*VCha(Cha4,1)*VChaC(Cha5,1) - VCha(Cha4,2)*
     -                    VChaC(Cha5,2)
                     bup49 = 2*SW2*IndexDelta(Cha4,Cha6) - 
     -                    2*VCha(Cha4,1)*VChaC(Cha6,1) - VCha(Cha4,2)*
     -                    VChaC(Cha6,2)
                     bup50 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -                    CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                     bup51 = -(cI*(CA*ZHggs(HiggsType,1))) - cI*(SA*
     -                    ZHggs(HiggsType,2)) + CB*ZHggs(HiggsType,3)
                     bup52 = cI*(SB*ZHggs(HiggsType,3)) + 
     -                    SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                     bup53 = cI*(SA*ZHggs(HiggsType,1)) - cI*(CA*
     -                    ZHggs(HiggsType,2)) + SB*ZHggs(HiggsType,3)
                     bup54 = cI*(bup46*bup49*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(bup42*bup45*(bup52*
     -                    UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    bup47*bup48*(bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    bup43*bup44*(bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))
                     bup55 = cI*(bup43*bup46*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(bup42*bup47*(bup52*
     -                    UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    bup45*bup48*(bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    bup44*bup49*(bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))
                     bup56 = cI*(bup43*bup44*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(bup47*bup48*(bup52*
     -                    UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    bup42*bup45*(bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    bup46*bup49*(bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))

                     loopa1 = loopa1 + (Alfa*cI)/(4.D0*CW2*Pi*SW2)*
     -                    (Cval(cc00,iintB18(Cha4,Cha5,Cha6))*
     -                    (-(bup56*EL*MCha(Cha5))-bup54*EL*MCha(Cha6)))/
     -                    (SW*Sqrt(2.D0)) + (Alfa*cI)/(8.D0*CW2*Pi*SW2)*
     -                    (Bval(bb0,iintB19(Cha5,Cha6))*(-(bup55*EL*
     -                    MCha(Cha4)) + bup56*EL*MCha(Cha5) + 
     -                    bup54*EL*MCha(Cha6)))/(SW*Sqrt(2.D0)) + 
     -                    Cval(cc2,iintB18(Cha4,Cha5,Cha6))*
     -                    ((Alfa*cI)/(8.D0*CW2*Pi*SW2)*
     -                    (q22*bup56*EL*MCha(Cha5))/(SW*Sqrt(2.D0)) + 
     -                    (Alfa*cI)/(16.D0*CW2*Pi*SW2)*(EL*MCha(Cha6)*
     -                    (cI*((q12+q22-q32)*bup46*bup49*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*((q12+q22-q32)*bup42*bup45*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    (q12+q22-q32)*bup47*bup48*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    (q12+q22-q32)*bup43*bup44*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/(SW*
     -                    Sqrt(2.D0)) - (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                    (EL*MCha(Cha4)*(cI*((q12+3*q22-q32)*bup43*
     -                    bup46*(bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*((q12+3*q22-q32)*bup42*bup47*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    (q12+3*q22-q32)*bup45*bup48*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    (q12+3*q22-q32)*bup44*bup49*(bup53*
     -                    UCha(Cha5,2)*VCha(Cha6,1) +bup51*UCha(Cha5,1)*
     -                    VCha(Cha6,2))))/(SW*Sqrt(2.D0)))

                     loopa1 = loopa1 +Cval(cc0,iintB18(Cha4,Cha5,Cha6))*
     -                    (-((Alfa*cI)/(8.D0*CW2*Pi*SW2)*
     -                    (bup55*EL*MCha(Cha4)*MCha2(Cha4))/
     -                    (SW*Sqrt(2.D0))) +(Alfa*cI)/(8.D0*CW2*Pi*SW2)*
     -                    (bup56*EL*MCha(Cha5)*MCha2(Cha4))/
     -                    (SW*Sqrt(2.D0))+(Alfa*cI)/(8.D0*CW2*Pi*SW2)*
     -                    (bup54*EL*MCha(Cha6)*MCha2(Cha4))/
     -                    (SW*Sqrt(2.D0))-(Alfa*cI)/(8.D0*CW2*Pi*SW2)*
     -                    (EL*MCha(Cha4)*MCha(Cha5)*MCha(Cha6)*
     -                    (cI*(bup44*bup49*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(bup45*bup48*(bup52*
     -                    UChaC(Cha5,2)*VChaC(Cha6,1) - bup50*
     -                    UChaC(Cha5,1)*VChaC(Cha6,2))) + bup42*bup47*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) +bup43*bup46*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + bup51*
     -                    UCha(Cha5,1)*VCha(Cha6,2))))/(SW*Sqrt(2.D0)) - 
     -                    (Alfa*cI)/(16.D0*CW2*Pi*SW2)*(EL*MCha(Cha4)*
     -                    (cI*((q12+q22-q32)*bup43*bup46*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*((q12+q22-q32)*bup42*bup47*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    (q12+q22-q32)*bup45*bup48*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    (q12+q22-q32)*bup44*bup49*(bup53*UCha(Cha5,2)*
     -                    VCha(Cha6,1) + bup51*UCha(Cha5,1)*
     -                    VCha(Cha6,2))))/(SW*Sqrt(2.D0))) + 
     -                    Cval(cc1,iintB18(Cha4,Cha5,Cha6))*((Alfa*cI)/
     -                    (8.D0*CW2*Pi*SW2)*(q12*bup54*EL*MCha(Cha6))/
     -                    (SW*Sqrt(2.D0)) - (3*Alfa*cI)/(16.D0*CW2*Pi*
     -                    SW2)*(EL*MCha(Cha4)*
     -                    (cI*((q12+1/3.D0*q22-1/3.D0*q32)*bup43*bup46*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*((q12+1/3.D0*q22-1/3.D0*q32)*bup42*bup47*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    (q12+1/3.D0*q22-1/3.D0*q32)*bup45*bup48*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    (q12+1/3.D0*q22-1/3.D0*q32)*bup44*bup49*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0)) +(Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                    (EL*MCha(Cha5)*(cI*((q12+q22-q32)*bup43*bup44*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*((q12+q22-q32)*bup47*bup48*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    (q12+q22-q32)*bup42*bup45*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    (q12+q22-q32)*bup46*bup49*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0)))


                     loopa2 = loopa2 + (Alfa*cI)/(4.D0*CW2*Pi*SW2)*
     -                    (Cval(cc12,iintB18(Cha4,Cha5,Cha6))*
     -                    (bup56*EL*MCha(Cha5) + bup54*EL*MCha(Cha6)))/
     -                    (SW*Sqrt(2.D0)) + 
     -                    Cval(cc2,iintB18(Cha4,Cha5,Cha6))*
     -                    ((Alfa*cI)/(16.D0*CW2*Pi*SW2)*(EL*MCha(Cha6)*
     -                    (cI*(2*bup46*bup49*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(2*bup42*bup45*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 2*bup47*
     -                    bup48*(bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 2*bup43*
     -                    bup44*(bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0)) -(Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                    (EL*MCha(Cha4)*(cI*(2*bup43*bup46*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*(2*bup42*bup47*(bup52*UChaC(Cha5,2)*
     -                    VChaC(Cha6,1) - bup50*UChaC(Cha5,1)*
     -                    VChaC(Cha6,2))) + 2*bup45*bup48*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 2*bup44*
     -                    bup49*(bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0)))
                     
                     loopa2 = loopa2+ Cval(cc0,iintB18(Cha4,Cha5,Cha6))*
     -                    (-(Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                    (EL*MCha(Cha4)*(cI*(2*bup43*bup46*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*(2*bup42*bup47*(bup52*UChaC(Cha5,2)*
     -                    VChaC(Cha6,1) - bup50*UChaC(Cha5,1)*
     -                    VChaC(Cha6,2))) + 2*bup45*bup48*(bup53*
     -                    UCha(Cha6,2)*VCha(Cha5,1) +bup51*UCha(Cha6,1)*
     -                    VCha(Cha5,2)) + 2*bup44*bup49*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0))) + 
     -                    Cval(cc1,iintB18(Cha4,Cha5,Cha6))*
     -                    (-(3*Alfa*cI)/(16.D0*CW2*Pi*SW2)*(EL*
     -                    MCha(Cha4)*(cI*((2/3D0)*bup43*bup46*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*((2/3D0)*bup42*bup47*(bup52*UChaC(Cha5,2)*
     -                    VChaC(Cha6,1) - bup50*UChaC(Cha5,1)*
     -                    VChaC(Cha6,2))) + (2/3D0)*bup45*bup48*(bup53*
     -                    UCha(Cha6,2)*VCha(Cha5,1) +bup51*UCha(Cha6,1)*
     -                    VCha(Cha5,2)) + (2/3D0)*bup44*bup49*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0)) +(Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                    (EL*MCha(Cha5)*(cI*(2*bup43*bup44*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*(2*bup47*bup48*(bup52*UChaC(Cha5,2)*
     -                    VChaC(Cha6,1) - bup50*UChaC(Cha5,1)*
     -                    VChaC(Cha6,2))) + 2*bup42*bup45*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    2*bup46*bup49*(bup53*UCha(Cha5,2)*
     -                    VCha(Cha6,1)+ bup51*UCha(Cha5,1)*
     -                    VCha(Cha6,2))))/(SW*Sqrt(2.D0)))


                     loopa3 = loopa3 +Cval(cc2,iintB18(Cha4,Cha5,Cha6))*
     -                    ((Alfa*cI)/(16.D0*CW2*Pi*SW2)*(EL*MCha(Cha6)*
     -                    (cI*(2*bup46*bup49*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(2*bup42*bup45*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + (-2)*
     -                    bup47*bup48*(bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + (-2)*bup43*
     -                    bup44*(bup53*UCha(Cha5,2)*VCha(Cha6,1) +bup51*
     -                    UCha(Cha5,1)*VCha(Cha6,2))))/(SW*Sqrt(2.D0)) - 
     -                    (Alfa*cI)/(16.D0*CW2*Pi*SW2)*(EL*MCha(Cha4)*
     -                    (cI*(2*bup43*bup46*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(2*bup42*bup47*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    (-2)*bup45*bup48*(bup53*UCha(Cha6,2)*
     -                    VCha(Cha5,1) + bup51*UCha(Cha6,1)*
     -                    VCha(Cha5,2)) + (-2)*bup44*bup49*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0)))
                     
                     loopa3 = loopa3 +Cval(cc0,iintB18(Cha4,Cha5,Cha6))*
     -                    (-(Alfa*cI)/(16.D0*CW2*Pi*SW2)*(EL*MCha(Cha4)*
     -                    (cI*(2*bup43*bup46*(bup52*UChaC(Cha6,2)*
     -                    VChaC(Cha5,1) - bup50*UChaC(Cha6,1)*
     -                    VChaC(Cha5,2))) + cI*(2*bup42*bup47*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + (-2)*
     -                    bup45*bup48*(bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + (-2)*bup44*
     -                    bup49*(bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0))) + 
     -                    Cval(cc1,iintB18(Cha4,Cha5,Cha6))*(-(3*Alfa*
     -                    cI)/(16.D0*CW2*Pi*SW2)*(EL*MCha(Cha4)*
     -                    (cI*((2/3D0)*bup43*bup46*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*((2/3D0)*bup42*bup47*
     -                    (bup52*UChaC(Cha5,2)*VChaC(Cha6,1) - 
     -                    bup50*UChaC(Cha5,1)*VChaC(Cha6,2))) + 
     -                    (-2/3D0)*bup45*bup48*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + 
     -                    (-2/3D0)*bup44*bup49*
     -                    (bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/(SW*
     -                    Sqrt(2.D0)) + (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                    (EL*MCha(Cha5)*(cI*(2*bup43*bup44*
     -                    (bup52*UChaC(Cha6,2)*VChaC(Cha5,1) - 
     -                    bup50*UChaC(Cha6,1)*VChaC(Cha5,2))) + 
     -                    cI*(2*bup47*bup48*(bup52*UChaC(Cha5,2)*
     -                    VChaC(Cha6,1) - bup50*UChaC(Cha5,1)*
     -                    VChaC(Cha6,2))) + (-2)*bup42*bup45*
     -                    (bup53*UCha(Cha6,2)*VCha(Cha5,1) + 
     -                    bup51*UCha(Cha6,1)*VCha(Cha5,2)) + (-2)*bup46*
     -                    bup49*(bup53*UCha(Cha5,2)*VCha(Cha6,1) + 
     -                    bup51*UCha(Cha5,1)*VCha(Cha6,2))))/
     -                    (SW*Sqrt(2.D0)))
                  enddo
               enddo
            enddo

            do Neu6 = 1,4
               do Neu5 = 1,4
                  do Neu4 = 1,4
                     bup4 = -(cI*(CB*ZHggs(HiggsType,3))) + CA*
     -                    ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
                     bup5 = cI*(CB*ZHggs(HiggsType,3)) + 
     -                    CA*ZHggs(HiggsType,1) +SA*ZHggs(HiggsType,2)
                     bup6 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -                    SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                     bup7 = cI*(SB*ZHggs(HiggsType,3)) + 
     -                    SA*ZHggs(HiggsType,1) -CA*ZHggs(HiggsType,2)
                     bup8 = -(CW*ZNeu(Neu5,2)) + SW*ZNeu(Neu5,1)
                     bup9 = -(CW*ZNeu(Neu6,2)) + SW*ZNeu(Neu6,1)
                     bup10 = ZNeu(Neu5,3)*ZNeuC(Neu4,3) - ZNeu(Neu5,4)*
     -                    ZNeuC(Neu4,4)
                     bup11 = ZNeu(Neu6,3)*ZNeuC(Neu4,3) - ZNeu(Neu6,4)*
     -                    ZNeuC(Neu4,4)
                     bup12 = -(CW*ZNeuC(Neu5,2)) + SW*ZNeuC(Neu5,1)
                     bup13 = ZNeu(Neu4,3)*ZNeuC(Neu5,3) - ZNeu(Neu4,4)*
     -                    ZNeuC(Neu5,4)
                     bup14 = -(CW*ZNeuC(Neu6,2)) + SW*ZNeuC(Neu6,1)
                     bup15 = ZNeu(Neu4,3)*ZNeuC(Neu6,3) - ZNeu(Neu4,4)*
     -                    ZNeuC(Neu6,4)
                     bup16 = bup6*bup9*ZNeu(Neu5,3) + 
     -                    bup5*bup9*ZNeu(Neu5,4) + 
     -                    bup8*(bup6*ZNeu(Neu6,3) + bup5*ZNeu(Neu6,4))
                     bup17 = bup14*bup7*ZNeuC(Neu5,3) + 
     -                    bup14*bup4*ZNeuC(Neu5,4) + 
     -                    bup12*(bup7*ZNeuC(Neu6,3) + bup4*ZNeuC(Neu6,4))

                     loopa1 = loopa1 + 1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc2,iintB8(Neu4,Neu5,Neu6))*
     -                    (Alfa*((q12+3*q22-q32)*bup13*bup15*bup16*EL*
     -                    MNeu(Neu4)) + Alfa*((q12+3*q22-q32)*bup10*
     -                    bup11*bup17*EL*MNeu(Neu4)) + 
     -                    Alfa*(q22*(2*bup10*bup15*bup16*EL*MNeu(Neu5) + 
     -                    2*bup11*bup13*bup17*EL*MNeu(Neu5))) + 
     -                    Alfa*((q12+q22-q32)*bup11*bup13*bup16*EL*
     -                    MNeu(Neu6)) + Alfa*((q12+q22-q32)*bup10*bup15*
     -                    bup17*EL*MNeu(Neu6))))/SW+1/(32.D0*CW*CW2*Pi*
     -                    SW2)*(Cval(cc1,iintB8(Neu4,Neu5,Neu6))*
     -                    (3*Alfa*((q12+1/3.D0*q22-1/3.D0*q32)*bup13*
     -                    bup15*bup16*EL*MNeu(Neu4)) + 
     -                    3*Alfa*((q12+1/3.D0*q22-1/3.D0*q32)*bup10*
     -                    bup11*bup17*EL*MNeu(Neu4)) + 
     -                    Alfa*((q12+q22-q32)*bup10*bup15*bup16*EL*
     -                    MNeu(Neu5)) + Alfa*((q12+q22-q32)*bup11*bup13*
     -                    bup17*EL*MNeu(Neu5)) + 
     -                    Alfa*(q12*(2*bup11*bup13*bup16*EL*MNeu(Neu6) + 
     -                    2*bup10*bup15*bup17*EL*MNeu(Neu6)))))/SW + 
     -                    Alfa/(8.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc00,iintB8(Neu4,Neu5,Neu6))*
     -                    (bup11*bup13*(-(bup17*EL*MNeu(Neu5)) - 
     -                    bup16*EL*MNeu(Neu6)) +bup10*bup15*(-(bup16*EL*
     -                    MNeu(Neu5)) - bup17*EL*MNeu(Neu6))))/SW + 
     -                    Alfa/(16.D0*CW*CW2*Pi*SW2)*
     -                    (Bval(bb0,iintB9(Neu5,Neu6))*
     -                    (bup13*bup15*bup16*EL*MNeu(Neu4) +bup11*bup13*
     -                    (bup17*EL*MNeu(Neu5) + bup16*EL*MNeu(Neu6)) + 
     -                    bup10*(bup11*bup17*EL*MNeu(Neu4) + 
     -                    bup15*(bup16*EL*MNeu(Neu5) + 
     -                    bup17*EL*MNeu(Neu6)))))/SW

                     loopa1 = loopa1 + 1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc0,iintB8(Neu4,Neu5,Neu6))*
     -                    (Alfa*((q12+q22-q32)*bup13*bup15*bup16*EL*
     -                    MNeu(Neu4)) + Alfa*((q12+q22-q32)*bup10*bup11*
     -                    bup17*EL*MNeu(Neu4)) + Alfa*(bup13*bup15*
     -                    (2*bup17*EL*MNeu(Neu4)*MNeu(Neu5)*MNeu(Neu6) + 
     -                    2*bup16*EL*MNeu(Neu4)*MNeu2(Neu4)) + bup11*
     -                    bup13*(2*bup17*EL*MNeu(Neu5)*MNeu2(Neu4) + 
     -                    2*bup16*EL*MNeu(Neu6)*MNeu2(Neu4)) + bup10*
     -                    (bup11*(2*bup16*EL*MNeu(Neu4)*MNeu(Neu5)*
     -                    MNeu(Neu6) + 2*bup17*EL*MNeu(Neu4)*
     -                    MNeu2(Neu4)) + bup15*
     -                    (2*bup16*EL*MNeu(Neu5)*MNeu2(Neu4) + 
     -                    2*bup17*EL*MNeu(Neu6)*MNeu2(Neu4))))))/SW


                     loopa2 = loopa2 + 1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc2,iintB8(Neu4,Neu5,Neu6))*
     -                    (Alfa*(2*bup13*bup15*bup16*EL*MNeu(Neu4)) + 
     -                    Alfa*(2*bup10*bup11*bup17*EL*MNeu(Neu4)) + 
     -                    Alfa*(2*bup11*bup13*bup16*EL*MNeu(Neu6)) + 
     -                    Alfa*(2*bup10*bup15*bup17*EL*MNeu(Neu6))))/SW+
     -                    1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc1,iintB8(Neu4,Neu5,Neu6))*(3*Alfa*
     -                    ((2/3D0)*bup13*bup15*bup16*EL*MNeu(Neu4)) + 3*
     -                    Alfa*((2/3D0)*bup10*bup11*bup17*EL*
     -                    MNeu(Neu4)) + Alfa*(2*bup10*bup15*bup16*EL*
     -                    MNeu(Neu5)) + Alfa*(2*bup11*bup13*bup17*EL*
     -                    MNeu(Neu5))))/SW + Alfa/(8.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc12,iintB8(Neu4,Neu5,Neu6))*(bup11*
     -                    bup13*(bup17*EL*MNeu(Neu5) + bup16*EL*
     -                    MNeu(Neu6))+bup10*bup15*(bup16*EL*MNeu(Neu5) + 
     -                    bup17*EL*MNeu(Neu6))))/SW 

                     loopa2 = loopa2 + 1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc0,iintB8(Neu4,Neu5,Neu6))*
     -                    (Alfa*(2*bup13*bup15*bup16*EL*MNeu(Neu4)) + 
     -                    Alfa*(2*bup10*bup11*bup17*EL*MNeu(Neu4))))/SW


                     loopa3 = loopa3 + 1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc2,iintB8(Neu4,Neu5,Neu6))*
     -                    (Alfa*((-2)*bup13*bup15*bup16*EL*MNeu(Neu4)) + 
     -                    Alfa*(2*bup10*bup11*bup17*EL*MNeu(Neu4)) + 
     -                    Alfa*((-2)*bup11*bup13*bup16*EL*MNeu(Neu6)) + 
     -                    Alfa*(2*bup10*bup15*bup17*EL*MNeu(Neu6))))/SW+
     -                    1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc1,iintB8(Neu4,Neu5,Neu6))*(3*Alfa*
     -                    ((-2/3D0)*bup13*bup15*bup16*EL*MNeu(Neu4)) + 
     -                    3*Alfa*((2/3D0)*bup10*bup11*bup17*EL*
     -                    MNeu(Neu4)) + Alfa*((-2)*bup10*bup15*bup16*EL*
     -                    MNeu(Neu5)) + Alfa*(2*bup11*bup13*bup17*EL*
     -                    MNeu(Neu5))))/SW 
                     
                     loopa3 = loopa3 + 1/(32.D0*CW*CW2*Pi*SW2)*
     -                    (Cval(cc0,iintB8(Neu4,Neu5,Neu6))*
     -                    (Alfa*((-2)*bup13*bup15*bup16*EL*MNeu(Neu4)) + 
     -                    Alfa*(2*bup10*bup11*bup17*EL*MNeu(Neu4))))/SW
                  enddo
               enddo
            enddo

** End of neutralino / chargino contributions
         end if
         
      end if
** end of 'bosonic' contribution


** The eps part of the loop - the second part converts between FC notation
** and standard notation (in FC Eps = -i*levi-civita)
      loopa3 = loopa3*(0.D0,-1.D0)

      if (debug .ne. 0) then
         write(*,*)'ZZh loopa1=', loopa1
         write(*,*)'ZZh loopa2=', loopa2
         write(*,*)'ZZh loopa3 =', loopa3
      end if


** Now we want the form factors themselves, for the light higgs
      a2vtxZ = -loopa2
         
      a1vtxZ = loopa1 + (loopa2*q1q2) + ct 
         
      a3vtxZ = loopa3
         
      if (debug .ne. 0) then
         write(*,*)"ZZh vtx a1 =", a1vtxZ
         write(*,*)"ZZh vtx a2 =", a2vtxZ
         write(*,*)"ZZh vtx a3 =", a3vtxZ
         write(*,*)"ZZh FERMION contn to vtx a1 =", a1vtxZ
         write(*,*)"ZZh FERMION contn to vtx a2 =", a2vtxZ
         write(*,*)"ZZh FERMION contn to vtx a3 =", a3vtxZ
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the photon-photon-
*** higgs vertex

      subroutine aa_vtxSUSY(a1vtxA,a2vtxA,a3vtxA,q12,q22,q32,q1q2,debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Dummy variables used when calculating the loop elements 
      double complex qup1, qup2
      double complex qup10, qup11, qup12, qup13, qup14, qup15, qup16
      double complex qup17, qup18, qup3, qup4, qup5, qup6, qup7
      double complex qup8, qup9
      double complex lup1
      double complex lup10, lup11, lup2, lup3, lup4, lup5, lup6
      double complex lup7, lup8, lup9
      double complex bup1, bup2, bup3
      double complex bup4, bup5

** The arguments of the integrals (suffix L -> leptons, B -> bosons)
      integer*8 iintQ1(3), iintQ2(3), iintQ3(3), iintQ4(3), iintQ5(3,2)
      integer*8 iintQ6(3,2), iintQ7(3,2), iintQ8(3,2)
      integer*8 iintL1(3), iintL2(3), iintL3(3,2), iintL4(3,2)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6, iintB7(2)
      integer*8 iintB8(2)

** Do loop variables
      integer Hig4, Hig5, Gen4, g4ini, sfe4, sfe5, Cha4, Neu4, Neu5

** These are the contributions from the counter-term (g-mu-nu), the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex loopa1, loopa2, loopa3


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita) 
** for photon-photon-higgs.  The first entry in each array is the light higgs, 
** the second entry is the heavy higgs and the third the cp-odd higgs
      double complex a1vtxA, a2vtxA, a3vtxA


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions
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


** Debugging, higgs and contributions parameter
      double complex UHggs(3,3), ZHggs(3,3)

** Identity matrix and imaginary i
      double precision IndexDelta(4,4)
      double precision Alfa, Alfa2

      integer debug


** Model parameters
#include "mssm.inc"      
#include "global.inc"

** This is needed for looptools
#include "looptools.h"


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions, 3 for SM-type particles
** plus sfermions and 4, 5, 6 for all particles
      if (sector .eq. 1) then
         g4ini = 3
      else 
         g4ini = 1
      end if


** UHggs and ZHggs: set to 1 for debugging
      do Hig4 = 1, 3
         do Hig5 = 1, 3
            if (Hig4 .eq. Hig5) then
               UHggs(Hig4,Hig5) = (1D0,0D0)
*               ZHggs(Hig4,Hig5) = (1D0,0D0)
            else
               UHggs(Hig4,Hig5) = (0D0,0D0)
*               ZHggs(Hig4,Hig5) = (0D0,0D0)
            end if
         end do
      end do


      if (propLEVEL .ne. 2) then
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               if (Hig4 .eq. Hig5) then
                  ZHggs(Hig4,Hig5) = (1D0,0D0)
               else
                  ZHggs(Hig4,Hig5) = (0D0,0D0)
               end if
            end do
         end do
      else
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               ZHggs(Hig4,Hig5) = Hgs(Hig4,Hig5)
            end do
         end do
      end if


** Clearing the values of the formfactors
      a1vtxA = (0D0,0D0)
      a2vtxA = (0D0,0D0)
      a3vtxA = (0D0,0D0)


** Setting identity matrix
      do Neu4 = 1, 4
         do Neu5 = 1, 4
            if (Neu4 .eq. Neu5) then
               IndexDelta(Neu4,Neu5) = 1D0
            else
               IndexDelta(Neu4,Neu5) = 0D0
            end if
         end do
      end do

      Alfa = AlfaQED
      Alfa2 = Alfa**2


** Now for the arguments of the integrals.  

** First the quarks and squarks
      do Gen4 = 1,3
         iintQ1(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ2(Gen4) = Cget(q22,q32,q12,Mf2(4,Gen4),
     -        Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ3(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
         iintQ4(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),
     -        Mf2(3,Gen4),Mf2(3,Gen4))
      enddo

      do Sfe4 = 1,2
         do Gen4 = 1,3
            iintQ5(Gen4,Sfe4) = 
     -           Bget(q32,MSf2(Sfe4,4,Gen4),MSf2(Sfe4,4,Gen4))
            iintQ6(Gen4,Sfe4) = 
     -           Cget(q22,q32,q12,MSf2(Sfe4,4,Gen4),MSf2(Sfe4,4,Gen4),
     -           MSf2(Sfe4,4,Gen4))
            iintQ7(Gen4,Sfe4) = 
     -           Bget(q32,MSf2(Sfe4,3,Gen4),MSf2(Sfe4,3,Gen4))
            iintQ8(Gen4,Sfe4) = 
     -           Cget(q22,q32,q12,MSf2(Sfe4,3,Gen4),MSf2(Sfe4,3,Gen4),
     -           MSf2(Sfe4,3,Gen4))
         enddo
      enddo

** The (s)leptonic integrals
      if (sector .ne. 1) then
         do Gen4 = 1,3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL2(Gen4) = Cget(q22,q32,q12,Mf2(2,Gen4),
     -           Mf2(2,Gen4),Mf2(2,Gen4))
         enddo
         do Sfe4 = 1,2
            do Gen4 = 1,3
               iintL3(Gen4,Sfe4) = 
     -              Bget(q32,MSf2(Sfe4,2,Gen4),MSf2(Sfe4,2,Gen4))
               iintL4(Gen4,Sfe4) = 
     -              Cget(q22,q32,q12,MSf2(Sfe4,2,Gen4),MSf2(Sfe4,2,Gen4),
     -              MSf2(Sfe4,2,Gen4))
            enddo
         enddo
      end if

** Bosonic integrals and their superpartners
      if (sector .ge. 3) then
         iintB1 = Bget(q22,MW2,MW2)
         iintB2 = Bget(q12,MW2,MW2)
         iintB3 = Bget(q32,MHp2,MHp2)
         iintB4 = Cget(q22,q32,q12,MHp2,MHp2,MHp2)
         iintB5 = Cget(q22,q32,q12,MW2,MW2,MW2)
         iintB6 = Bget(q32,MW2,MW2)

         if (sector .ge. 4) then
            do Cha4 = 1,2
               iintB7(Cha4) = Bget(q32,MCha2(Cha4),MCha2(Cha4))
               iintB8(Cha4) = Cget(q22,q32,q12,MCha2(Cha4),
     -              MCha2(Cha4),MCha2(Cha4))
            enddo
         end if
      end if


** Clearing the values of each contribution
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0) 

** Now the insanely long and hideously complicated form factors.  In bits.  
** Now the loop

** Quark and squark contributions
      qup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
      qup2 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)

      do Gen4 = 1,3
         loopa1 = loopa1 + Alfa/(MW*Pi)*(EL*(((qup2*(
     -        (4/3.D0*Bval(bb0,iintQ3(Gen4)) - 
     -        16/3.D0*Cval(cc00,iintQ4(Gen4)))) - 
     -        1/3.D0*(Cval(cc0,iintQ4(Gen4))*((q12+q22-q32)*
     -        (-(cI*(SB*ZHggs(HiggsType,3))) + 
     -        TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)))+
     -        (q12+q22-q32)*(cI*(SB*ZHggs(HiggsType,3)) + 
     -        TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))))
     -        )/TB)*Mf2(3,Gen4))/SB + ((qup1*(
     -        (-(1/3.D0*Bval(bb0,iintQ1(Gen4))) + 
     -        4/3.D0*Cval(cc00,iintQ2(Gen4)))) + 
     -        1/12.D0*(Cval(cc0,iintQ2(Gen4))*((q12+q22-q32)*
     -        (-(cI*(SB*ZHggs(HiggsType,3))) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)) + (q12+q22-q32)*
     -        (cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)))))*Mf2(4,Gen4))/CB))/SW

         loopa2 = loopa2 + Alfa/(MW*Pi)*(EL*(((qup2*(16/3.D0*
     -        (Cval(cc12,iintQ4(Gen4)))) - 1/3.D0*
     -        (Cval(cc0,iintQ4(Gen4))*(2*(-(cI*(SB*ZHggs(HiggsType,3)))+ 
     -        TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)))+2*
     -        (cI*(SB*ZHggs(HiggsType,3)) + TB*(CA*ZHggs(HiggsType,1) + 
     -        SA*ZHggs(HiggsType,2)))))/TB)*Mf2(3,Gen4))/SB + ((qup1*
     -        (-4/3.D0*(Cval(cc12,iintQ2(Gen4)))) + 1/12.D0*
     -        (Cval(cc0,iintQ2(Gen4))*(2*(-(cI*(SB*ZHggs(HiggsType,3)))+ 
     -        SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 2*
     -        (cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)))))*Mf2(4,Gen4))/CB))/SW

         loopa3 = loopa3 + Alfa/(MW*Pi)*(EL*(((-1/3.D0*
     -        (Cval(cc0,iintQ4(Gen4))*((-2)*(-(cI*(SB*
     -        ZHggs(HiggsType,3))) + TB*(CA*ZHggs(HiggsType,1) + SA*
     -        ZHggs(HiggsType,2)))+2*(cI*(SB*ZHggs(HiggsType,3)) + 
     -        TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))))
     -        )/TB)*Mf2(3,Gen4))/SB + ((1/12.D0*(Cval(cc0,iintQ2(Gen4))*
     -        (2*(-(cI*(SB*ZHggs(HiggsType,3))) +SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)) + (-2)*
     -        (cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)))))*Mf2(4,Gen4))/CB))/SW
      enddo

      qup3 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
      qup4 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
      qup5 = -(SAB*ZHggs(HiggsType,1)) + CAB*ZHggs(HiggsType,2)
      qup6 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
      qup7 = -(cI*(CB*ZHggs(HiggsType,3))) + CA*ZHggs(HiggsType,1) + 
     -     SA*ZHggs(HiggsType,2)
      qup8 = cI*(CB*ZHggs(HiggsType,3)) + CA*ZHggs(HiggsType,1) + 
     -     SA*ZHggs(HiggsType,2)
      qup9 = -(cI*(SB*ZHggs(HiggsType,3))) + SA*ZHggs(HiggsType,1) - 
     -     CA*ZHggs(HiggsType,2)
      qup10 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -     TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))
      qup11 = cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -     CA*ZHggs(HiggsType,2)
      qup12 = cI*(SB*ZHggs(HiggsType,3)) + 
     -     TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))

      do Sfe4 = 1,2
         do Gen4 = 1,3
            qup13 = MW*MZ*(-3 + 4*SW2)*(qup4*SB)+6*CW*(qup6*Mf2(3,Gen4))
            qup14 = MW*MZ*(-3 + 2*SW2)*(CB*qup4)+6*CW*(qup3*Mf2(4,Gen4))
            qup15 = 3*CW*((qup9*MUE*TB + qup10*AfC(3,Gen4,Gen4))*
     -           Mf(3,Gen4)*USfC(Sfe4,2,3,Gen4)) + 
     -           qup13*TB*USfC(Sfe4,1,3,Gen4)
            qup16 = 3*CW*((qup11*MUEC*TB + qup12*Af(3,Gen4,Gen4))*
     -           Mf(3,Gen4)*USfC(Sfe4,1,3,Gen4)) + 
     -           2*TB*(2*MW*MZ*SW2*(qup5*SB) + 3*CW*(qup6*Mf2(3,Gen4)))*
     -           USfC(Sfe4,2,3,Gen4)
            qup17 = 3*CW*((qup8*MUE + qup11*AfC(4,Gen4,Gen4))*
     -           Mf(4,Gen4)*USfC(Sfe4,2,4,Gen4)) + 
     -           qup14*USfC(Sfe4,1,4,Gen4)
            qup18 = 3*CW*((qup7*MUEC + qup9*Af(4,Gen4,Gen4))*
     -           Mf(4,Gen4)*USfC(Sfe4,1,4,Gen4)) + 
     -           2*(MW*MZ*SW2*(CB*qup5) + 3*CW*(qup3*Mf2(4,Gen4)))*
     -           USfC(Sfe4,2,4,Gen4)

            loopa1 = loopa1 + Alfa/(CW*MW*Pi)*(EL*
     -           ((((-(1/9.D0*Bval(bb0,iintQ7(Gen4,Sfe4))) + 
     -           4/9.D0*Cval(cc00,iintQ8(Gen4,Sfe4))))*
     -           (qup15*USf(Sfe4,1,3,Gen4) + 
     -           qup16*USf(Sfe4,2,3,Gen4)))/(SB*TB) + 
     -           (((1/36.D0*Bval(bb0,iintQ5(Gen4,Sfe4)) - 
     -           1/9.D0*Cval(cc00,iintQ6(Gen4,Sfe4))))*
     -           (qup17*USf(Sfe4,1,4,Gen4) + 
     -           qup18*USf(Sfe4,2,4,Gen4)))/CB))/SW

            loopa2 = loopa2 + Alfa/(CW*MW*Pi)*(EL*
     -           (((-4/9.D0*(Cval(cc12,iintQ8(Gen4,Sfe4))))*
     -           (qup15*USf(Sfe4,1,3,Gen4) + qup16*USf(Sfe4,2,3,Gen4)))/
     -           (SB*TB) + ((1/9.D0*(Cval(cc12,iintQ6(Gen4,Sfe4))))*
     -           (qup17*USf(Sfe4,1,4,Gen4) + 
     -           qup18*USf(Sfe4,2,4,Gen4)))/CB))/SW
         enddo
      enddo


** The lepton and slepton contributions
      if (sector .ne. 1) then
         lup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
         do Gen4 = 1,3
            loopa1 = loopa1 + Alfa/(MW*Pi)*(EL*(lup1*
     -           ((-Bval(bb0,iintL1(Gen4)) +4*Cval(cc00,iintL2(Gen4))))+ 
     -           1/4.D0*(Cval(cc0,iintL2(Gen4))*((q12+q22-q32)*
     -           (-(cI*(SB*ZHggs(HiggsType,3))) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)) + (q12+q22-q32)*
     -           (cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)))))*Mf2(2,Gen4))/(CB*SW)

            loopa2 = loopa2 + Alfa/(MW*Pi)*(EL*(lup1*(-4*
     -           Cval(cc12,iintL2(Gen4))) + 1/4.D0*
     -           (Cval(cc0,iintL2(Gen4))*(2*(-(cI*(SB*
     -           ZHggs(HiggsType,3))) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)) + 2*
     -           (cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)))))*Mf2(2,Gen4))/(CB*SW)

            loopa3 = loopa3 + Alfa/(MW*Pi)*(EL*(1/4.D0*
     -           (Cval(cc0,iintL2(Gen4))*(2*(-(cI*(SB*
     -           ZHggs(HiggsType,3))) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)) - 2*
     -           (cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)))))*Mf2(2,Gen4))/(CB*SW)
         enddo


         lup2 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
         lup3 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
         lup4 = -(SAB*ZHggs(HiggsType,1)) + CAB*ZHggs(HiggsType,2)
         lup5 = -(cI*(CB*ZHggs(HiggsType,3))) + CA*ZHggs(HiggsType,1) + 
     -        SA*ZHggs(HiggsType,2)
         lup6 = cI*(CB*ZHggs(HiggsType,3)) + CA*ZHggs(HiggsType,1) + SA*
     -        ZHggs(HiggsType,2)
         lup7 = -(cI*(SB*ZHggs(HiggsType,3))) + SA*ZHggs(HiggsType,1) - 
     -        CA*ZHggs(HiggsType,2)
         lup8 = cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - CA*
     -        ZHggs(HiggsType,2)

         do Sfe4 = 1,2
            do Gen4 = 1,3
               lup9 = MW*MZ*(-1 + 2*SW2)*(CB*lup3) + 2*CW*(lup2*
     -              Mf2(2,Gen4))
               lup10 = CW*((lup6*MUE +lup8*AfC(2,Gen4,Gen4))*Mf(2,Gen4)*
     -              USfC(Sfe4,2,2,Gen4)) + lup9*USfC(Sfe4,1,2,Gen4)
               lup11 = CW*((lup5*MUEC + lup7*Af(2,Gen4,Gen4))*
     -              Mf(2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 2*(MW*MZ*SW2*(CB*
     -              lup4) + CW*(lup2*Mf2(2,Gen4)))*USfC(Sfe4,2,2,Gen4)

               loopa1 = loopa1 + Alfa/(CW*MW*Pi)*(EL*
     -              ((-Cval(cc00,iintL4(Gen4,Sfe4)) + 1/4.D0*
     -              Bval(bb0,iintL3(Gen4,Sfe4))))*(lup10*
     -              USf(Sfe4,1,2,Gen4) + lup11*USf(Sfe4,2,2,Gen4)))/
     -              (CB*SW)

               loopa2 = loopa2 + Alfa/(CW*MW*Pi)*(EL*
     -              (Cval(cc12,iintL4(Gen4,Sfe4)))*(lup10*
     -              USf(Sfe4,1,2,Gen4) + lup11*USf(Sfe4,2,2,Gen4)))/
     -              (CB*SW)
            enddo
         enddo
      end if


** fermion formfactors 
      a2vtxAHfm = -loopa2         
      a1vtxAHfm = loopa1 + (loopa2*q1q2)          
      a3vtxAHfm = loopa3


** And the SM bosons and their superpartners ...
      if (sector .ge. 3) then
         bup1 = SBA*ZHggs(HiggsType,1) + CBA*ZHggs(HiggsType,2)
         bup2 = 2*CW2*(SBA*ZHggs(HiggsType,1) +CBA*ZHggs(HiggsType,2)) + 
     -        C2B*(SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2))
         bup3 = 9*CW2*(SBA*ZHggs(HiggsType,1) +CBA*ZHggs(HiggsType,2)) + 
     -        C2B*(-2*SAB*ZHggs(HiggsType,1) + 2*CAB*ZHggs(HiggsType,2))

         loopa1 = loopa1 + (Alfa*MW)/Pi*(EL*(1/4.D0*(bup1*(((q12-q22-
     -        q32) + (-q32 + 3*q22 + q12) - (q12-q22-q32))*
     -        Cval(cc1,iintB5) + ((q12-q22+q32) + 3*(q12 - 1/3.D0*q32 + 
     -        1/3.D0*q22) - (q12-q22+q32))*Cval(cc2,iintB5))) + 
     -        1/CW2*(1/2.D0*(Cval(cc00,iintB5)*(bup3 + 3*CW2*bup1)) + 
     -        1/4.D0*(Cval(cc0,iintB5)*(MW2*(bup2) + CW2*(bup1*
     -        (5*(q12+q22-q32) + (q12-q32) + (q22-q32) - 2*MW2))))) + 
     -        1*(1/CW2*(bup2*Cval(cc00,iintB4)) + 1/4.D0*(-(bup1*
     -        (Bval(bb0,iintB1) + Bval(bb0,iintB2))) - 1/CW2*(bup2*
     -        Bval(bb0,iintB3))) + Bval(bb0,iintB6)*(-(1/2.D0*bup1) - 
     -        1/(4.D0*CW2)*(2*CW2*(SBA*ZHggs(HiggsType,1) + CBA*
     -        ZHggs(HiggsType,2)) + C2B*(-(SAB*ZHggs(HiggsType,1)) + 
     -        CAB*ZHggs(HiggsType,2)))))))/SW

         loopa2 = loopa2 + (Alfa*MW)/Pi*(EL*(1/CW2*(-(bup2*
     -        Cval(cc12,iintB4)) + 1/2.D0*(Cval(cc12,iintB5)*
     -        (-bup3 - 3*CW2*bup1)) + 1/4.D0*(Cval(cc0,iintB5)*(CW2*
     -        (bup1*(16D0)))))))/SW
         
** chargino / neutralino contribution
         if (sector .ge. 4) then

            bup4 = cI*(CB*ZHggs(HiggsType,3)) + CA*ZHggs(HiggsType,1) + 
     -           SA*ZHggs(HiggsType,2)
            bup5 = -(cI*(SB*ZHggs(HiggsType,3)))+SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)

            do Cha4 = 1,2
               loopa1 = loopa1 + Alfa/Pi*(EL*MCha(Cha4)*(1/2.D0*
     -              (Cval(cc0,iintB8(Cha4))*((q12+q22-q32)*(bup5*
     -              UCha(Cha4,2)*VCha(Cha4,1) - bup4*UCha(Cha4,1)*
     -              VCha(Cha4,2)) + (q12+q22-q32)*(UChaC(Cha4,2)*
     -              VChaC(Cha4,1)*(cI*(SB*ZHggs(HiggsType,3)) + SA*
     -              ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) - 
     -              UChaC(Cha4,1)*VChaC(Cha4,2)*(-(cI*(CB*
     -              ZHggs(HiggsType,3))) + CA*ZHggs(HiggsType,1) + SA*
     -              ZHggs(HiggsType,2)))))/Sqrt(2.D0) + 
     -              ((-(Bval(bb0,iintB7(Cha4))/Sqrt(2.D0)) + 
     -              2*Cval(cc00,iintB8(Cha4))*Sqrt(2.D0)))*
     -              (cI*((SB*UChaC(Cha4,2)*VChaC(Cha4,1) + CB*
     -              UChaC(Cha4,1)*VChaC(Cha4,2))*ZHggs(HiggsType,3)) +
     -              bup5*UCha(Cha4,2)*VCha(Cha4,1) - bup4*UCha(Cha4,1)*
     -              VCha(Cha4,2) + CA*(-(UChaC(Cha4,1)*VChaC(Cha4,2)*
     -              ZHggs(HiggsType,1)) - UChaC(Cha4,2)*VChaC(Cha4,1)*
     -              ZHggs(HiggsType,2)) + SA*(UChaC(Cha4,2)*
     -              VChaC(Cha4,1)*ZHggs(HiggsType,1) - UChaC(Cha4,1)*
     -              VChaC(Cha4,2)*ZHggs(HiggsType,2)))))/SW

               loopa2 = loopa2 + Alfa/Pi*(EL*MCha(Cha4)*(1/2.D0*
     -              (Cval(cc0,iintB8(Cha4))*(2*(bup5*UCha(Cha4,2)*
     -              VCha(Cha4,1) - bup4*UCha(Cha4,1)*VCha(Cha4,2)) + 2*
     -              (UChaC(Cha4,2)*VChaC(Cha4,1)*(cI*(SB*
     -              ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - CA*
     -              ZHggs(HiggsType,2)) - UChaC(Cha4,1)*VChaC(Cha4,2)*
     -              (-(cI*(CB*ZHggs(HiggsType,3))) + CA*
     -              ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)))))/
     -              Sqrt(2.D0) +(-2*Cval(cc12,iintB8(Cha4))*Sqrt(2.D0))*
     -              (cI*((SB*UChaC(Cha4,2)*VChaC(Cha4,1) + CB*
     -              UChaC(Cha4,1)*VChaC(Cha4,2))*ZHggs(HiggsType,3)) +
     -              bup5*UCha(Cha4,2)*VCha(Cha4,1) - bup4*UCha(Cha4,1)*
     -              VCha(Cha4,2) + CA*(-(UChaC(Cha4,1)*VChaC(Cha4,2)*
     -              ZHggs(HiggsType,1)) - UChaC(Cha4,2)*VChaC(Cha4,1)*
     -              ZHggs(HiggsType,2)) + SA*(UChaC(Cha4,2)*
     -              VChaC(Cha4,1)*ZHggs(HiggsType,1) - UChaC(Cha4,1)*
     -              VChaC(Cha4,2)*ZHggs(HiggsType,2)))))/SW

               loopa3 = loopa3 + Alfa/Pi*(EL*MCha(Cha4)*(1/2.D0*
     -              (Cval(cc0,iintB8(Cha4))*((-2)*(bup5*UCha(Cha4,2)*
     -              VCha(Cha4,1) - bup4*UCha(Cha4,1)*VCha(Cha4,2)) + 2*
     -              (UChaC(Cha4,2)*VChaC(Cha4,1)*(cI*(SB*
     -              ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - CA*
     -              ZHggs(HiggsType,2)) - UChaC(Cha4,1)*VChaC(Cha4,2)*
     -              (-(cI*(CB*ZHggs(HiggsType,3))) + CA*
     -              ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)))))/
     -              Sqrt(2.D0)))/SW
            enddo
         end if
      end if
** end of 'bosonic' contribution


** The eps part of the loop - the second part converts between FC notation
** and standard notation (in FC Eps = -i*levi-civita)
      loopa3 = loopa3*(0.D0,-1.D0)

      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'photon-photon-h loopa1 =', loopa1
         write(*,*)'photon-photon-h loopa2 =', loopa2
         write(*,*)'photon-photon-h loopa3 =', loopa3
      end if


** Now we want the form factors themselves
      a2vtxA = -loopa2         
      a1vtxA = loopa1 + (loopa2*q1q2)          
      a3vtxA = loopa3
         
      if (debug .ne. 0) then
         write(*,*)"photon-photon-h vtx a1 =", a1vtxA
         write(*,*)"photon-photon-h vtx a2 =", a2vtxA
         write(*,*)"photon-photon-h vtx a3 =", a3vtxA
         write(*,*)'FERMION contn to a1vtxA = ', a1vtxAHfm
         write(*,*)'FERMION contn to a2vtxA = ', a2vtxAHfm
         write(*,*)'FERMION contn to a3vtxA = ', a3vtxAHfm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the photon-Z-higgs 
*** vertex

      subroutine az_vtxSUSY(a1vtxAZ,a2vtxAZ,a3vtxAZ,q12,q22,q32,q1q2,
     -     debug,fm)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2

** Dummy variables used when calculating the loop elements
      double complex qup1, qup2, qup3, qup4, qup5, qup6
      double complex qup10, qup11, qup12, qup13, qup14, qup15, qup16
      double complex qup17, qup18, qup19, qup20, qup21, qup22, qup23
      double complex qup24, qup25, qup26, qup27, qup28, qup29, qup30
      double complex qup7, qup8, qup9
      double complex lup1, lup2, lup3
      double complex lup10, lup11, lup12, lup13, lup14, lup15, lup16
      double complex lup17, lup4, lup5, lup6, lup7, lup8, lup9
      double complex bup1, bup2, bup3
      double complex bup10, bup11, bup12, bup13, bup14, bup4, bup5
      double complex bup6, bup7, bup8, bup9

** The arguments of the integrals (suffix L -> lepton, B -> boson)
      integer*8 iintQ1(3), iintQ2(3), iintQ3(3), iintQ4(3), iintQ5(3,2,2)
      integer*8 iintQ6(3,2,2), iintQ7(3,2,2), iintQ8(3,2,2)
      integer*8 iintL1(3), iintL2(3), iintL3(3,2,2), iintL4(3,2,2)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6
      integer*8 iintB7(2,2), iintB8(2,2)

** Do loop variables
      integer g4ini, Gen4, Sfe5, Sfe4, Cha4, Cha5, Hig4, Hig5, Neu4, Neu5

** These are the contributions from the counter-term (g-mu-nu), the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex loopa1, loopa2, loopa3, ct


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita) 
** for photon-photon-higgs.  The first entry in each array is the light higgs, 
** the second entry is the heavy higgs and the third the cp-odd higgs.
      double complex a1vtxAZ, a2vtxAZ, a3vtxAZ


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions.  'fm' = 1 for AZ and 2 for AZ2
      integer fm
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


** Higgs and propagator contributions parameter
      double complex UHggs(3,3), ZHggs(3,3)


** Identity matrix and imaginary i
      double precision IndexDelta(4,4)
      double precision Alfa, Alfa2

      integer debug

** Model parameters
#include "mssm.inc"
#include "global.inc"

** This is needed for looptools
#include "looptools.h"


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions, 3 for SM-type particles
** plus sfermions and 4, 5, 6 for all particles
      if (sector .eq. 1) then
         g4ini = 3
      else 
         g4ini = 1
      end if


** UHggs and ZHggs: set to 1 for debugging
      do Hig4 = 1, 3
         do Hig5 = 1, 3
            if (Hig4 .eq. Hig5) then
               UHggs(Hig4,Hig5) = (1D0,0D0)
*               ZHggs(Hig4,Hig5) = (1D0,0D0)
            else
               UHggs(Hig4,Hig5) = (0D0,0D0)
*               ZHggs(Hig4,Hig5) = (0D0,0D0)
            end if
         end do
      end do


      if (propLEVEL .ne. 2) then
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               if (Hig4 .eq. Hig5) then
                  ZHggs(Hig4,Hig5) = (1D0,0D0)
               else
                  ZHggs(Hig4,Hig5) = (0D0,0D0)
               end if
            end do
         end do
      else
         do Hig4 = 1, 3
            do Hig5 = 1, 3
               ZHggs(Hig4,Hig5) = Hgs(Hig4,Hig5)
            end do
         end do
      end if


** Clearing the values of the formfactors
      a1vtxAZ = (0D0,0D0)
      a2vtxAZ = (0D0,0D0)
      a3vtxAZ = (0D0,0D0)


** Setting identity matrix
      do Neu4 = 1, 4
         do Neu5 = 1, 4
            if (Neu4 .eq. Neu5) then
               IndexDelta(Neu4,Neu5) = 1D0
            else
               IndexDelta(Neu4,Neu5) = 0D0
            end if
         end do
      end do

      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Now for the arguments of the integrals.  

** First the quarks and squarks
      do Gen4 = 1,3
         iintQ1(Gen4) = Cget(q22,q32,q12,Mf2(4,Gen4),
     -        Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ2(Gen4) = Cget(q22,q32,q12,Mf2(3,Gen4),
     -        Mf2(3,Gen4),Mf2(3,Gen4))
         iintQ3(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
         iintQ4(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
      enddo
      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = 1,3
               iintQ5(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,4,Gen4),MSf2(Sfe5,4,Gen4))
               iintQ6(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,3,Gen4))
               iintQ7(Gen4,Sfe4,Sfe5) = 
     -              Cget(q12,q32,q22,MSf2(Sfe4,4,Gen4),MSf2(Sfe4,4,Gen4),
     -              MSf2(Sfe5,4,Gen4))
               iintQ8(Gen4,Sfe4,Sfe5) = 
     -              Cget(q12,q32,q22,MSf2(Sfe4,3,Gen4),MSf2(Sfe4,3,Gen4),
     -              MSf2(Sfe5,3,Gen4))
            enddo
         enddo
      enddo

** The (s)leptonic integrals
      if (sector .ne. 1) then
         do Gen4 = 1,3
            iintL1(Gen4) = Cget(q22,q32,q12,Mf2(2,Gen4),
     -           Mf2(2,Gen4),Mf2(2,Gen4))
            iintL2(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
         enddo
         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  iintL3(Gen4,Sfe4,Sfe5) = 
     -                 Bget(q32,MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4))
                  iintL4(Gen4,Sfe4,Sfe5) = 
     -                 Cget(q12,q32,q22,MSf2(Sfe4,2,Gen4),
     -                 MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4))
               enddo
            enddo
         enddo
      end if

** Bosonic integrals and their superpartners
      if (sector .ge. 3) then
         iintB1 = Bget(q12,MW2,MW2)
         iintB2 = Bget(q22,MW2,MW2)
         iintB3 = Bget(q32,MHp2,MHp2)
         iintB4 = Cget(q22,q32,q12,MHp2,MHp2,MHp2)
         iintB5 = Cget(q22,q32,q12,MW2,MW2,MW2)
         iintB6 = Bget(q32,MW2,MW2)

** chargino/neutralino contributions         
         if (sector .ge. 4) then
            do Cha5 = 1,2
               do Cha4 = 1,2
                  iintB7(Cha4,Cha5) = Bget(q32,MCha2(Cha4),MCha2(Cha5))
                  iintB8(Cha4,Cha5) = Cget(q12,q32,q22,
     -                 MCha2(Cha4),MCha2(Cha4),MCha2(Cha5))
               enddo
            enddo
         end if
      end if


** Clearing the values of each contribution
      loopa1 = (0D0,0D0)  
      loopa2 = (0D0,0D0) 
      loopa3 = (0D0,0D0) 
      ct = (0D0,0D0)

** Now the insanely long and hideously complicated form factors.  In bits. 

** The counterterm:
* h:         ct = MW/(2.D0*CW2)*(dZZA1*EL*SBA)/SW 
      ct = MW/(2.D0*CW2)*(dZZA1*EL*(SBA*ZHggs(HiggsType,1) +CBA*
     -     ZHggs(HiggsType,2)))/SW

      if (sector .eq. 4) then
         ct = MW/(2.D0*CW2)*(dZZA1full*EL*(SBA*ZHggs(HiggsType,1) +
     -        CBA*ZHggs(HiggsType,2)))/SW
* h0:             ct = MW/(2.D0*CW2)*(dZZA1full*EL*SBA)/SW 
      end if

* for the fermions
         ctFM = MW/(2.D0*CW2)*(dZZA1fm*EL*SBA)/SW 

      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'photon-Z-h vertex counterterm =', ct
            write(*,*)'FERMION contn to photon-Z-h vertex counterterm =', ct
      end if

** Now the loops:

** Quark and squark contributions
      qup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
      qup2 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
      qup3 = -(cI*(SB*ZHggs(HiggsType,3))) + SA*ZHggs(HiggsType,1) - 
     -     CA*ZHggs(HiggsType,2)
      qup4 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -     TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))
      qup5 = cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -     CA*ZHggs(HiggsType,2)
      qup6 = cI*(SB*ZHggs(HiggsType,3)) + 
     -     TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))

      do Gen4 = 1,3
         loopa1 = loopa1 +Alfa/(CW*MW*Pi*SW2)*(EL*((((-3 + 8*SW2)/24.D0*
     -        (((-q32 + q22 + q12)*qup4 + (-q32 + q22 + q12)*qup6)*
     -        Cval(cc0,iintQ2(Gen4)))/TB + (3 - 8*SW2)*(qup2*((1/6.D0*
     -        Bval(bb0,iintQ4(Gen4))-2/3.D0*Cval(cc00,iintQ2(Gen4))))))*
     -        Mf2(3,Gen4))/SB + (((3 - 4*SW2)/48.D0*(((-q32 + q22 +q12)*
     -        qup3 + (-q32 + q22 + q12)*qup5)*Cval(cc0,iintQ1(Gen4))) + 
     -        qup1*((-3 + 4*SW2)/12.D0*(Bval(bb0,iintQ3(Gen4))) + 
     -        (3 - 4*SW2)*(1/3.D0*(Cval(cc00,iintQ1(Gen4))))))*
     -        Mf2(4,Gen4))/CB))

         loopa2 = loopa2 + Alfa/(CW*MW*Pi*SW2)*(EL*((((-3+ 8*SW2)/24.D0*
     -        ((2*qup4 + 2*qup6)*Cval(cc0,iintQ2(Gen4)))/TB + 
     -        (3 - 8*SW2)*(qup2*(2/3.D0*(Cval(cc12,iintQ2(Gen4))))))*
     -        Mf2(3,Gen4))/SB + (((3 - 4*SW2)/48.D0*((2*qup3 + 2*qup5)*
     -        Cval(cc0,iintQ1(Gen4))) + qup1*((3 - 4*SW2)*(-1/3.D0*
     -        (Cval(cc12,iintQ1(Gen4))))))*Mf2(4,Gen4))/CB))

         loopa3 = loopa3 + Alfa/(CW*MW*Pi*SW2)*(EL*((((-3 +8*SW2)/24.D0*
     -        (((-2)*qup4 + 2*qup6)*Cval(cc0,iintQ2(Gen4)))/TB)*
     -        Mf2(3,Gen4))/SB))
      enddo

      qup7 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
      qup8 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
      qup9 = -(SAB*ZHggs(HiggsType,1)) + CAB*ZHggs(HiggsType,2)
      qup10 = CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
      qup11 = -(cI*(CB*ZHggs(HiggsType,3))) + 
     -     CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
      qup12 = cI*(CB*ZHggs(HiggsType,3)) + 
     -     CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
      qup13 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -     SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
      qup14 = -(cI*(SB*ZHggs(HiggsType,3))) + 
     -     TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))
      qup15 = cI*(SB*ZHggs(HiggsType,3)) + 
     -     SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
      qup16 = cI*(SB*ZHggs(HiggsType,3)) + 
     -     TB*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = 1,3
               qup17 = MW*MZ*(-3 + 4*SW2)*(qup8*SB) + 6*CW*
     -              (qup10*Mf2(3,Gen4))
               qup18 = MW*MZ*(-3 + 2*SW2)*(CB*qup8) + 6*CW*
     -              (qup7*Mf2(4,Gen4))
               qup19 = 4*SW2*(USf(Sfe5,1,3,Gen4)*USfC(Sfe4,1,3,Gen4) + 
     -              USf(Sfe5,2,3,Gen4)*USfC(Sfe4,2,3,Gen4)) - 
     -              3*USf(Sfe5,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)
               qup20 = 2*SW2*(USf(Sfe5,1,4,Gen4)*USfC(Sfe4,1,4,Gen4) + 
     -              USf(Sfe5,2,4,Gen4)*USfC(Sfe4,2,4,Gen4)) - 
     -              3*USf(Sfe5,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)
               qup21 = 4*SW2*(USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4) + 
     -              USf(Sfe4,2,3,Gen4)*USfC(Sfe5,2,3,Gen4)) - 
     -              3*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)
               qup22 = 2*SW2*(USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4) + 
     -              USf(Sfe4,2,4,Gen4)*USfC(Sfe5,2,4,Gen4)) - 
     -              3*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)
               qup23 = 3*CW*((qup13*MUE*TB + qup14*AfC(3,Gen4,Gen4))*
     -              Mf(3,Gen4)*USfC(Sfe4,2,3,Gen4)) + 
     -              qup17*TB*USfC(Sfe4,1,3,Gen4)
               qup24 = 3*CW*((qup15*MUEC*TB + qup16*Af(3,Gen4,Gen4))*
     -              Mf(3,Gen4)*USfC(Sfe4,1,3,Gen4)) + 2*TB*(2*MW*MZ*SW2*
     -              (qup9*SB) + 3*CW*(qup10*Mf2(3,Gen4)))*
     -              USfC(Sfe4,2,3,Gen4)
               qup25 = 3*CW*((qup12*MUE + qup15*AfC(4,Gen4,Gen4))*
     -              Mf(4,Gen4)*USfC(Sfe4,2,4,Gen4)) + 
     -              qup18*USfC(Sfe4,1,4,Gen4)
               qup26 = 3*CW*((qup11*MUEC + qup13*Af(4,Gen4,Gen4))*
     -              Mf(4,Gen4)*USfC(Sfe4,1,4,Gen4)) + 
     -              2*(MW*MZ*SW2*(CB*qup9) + 3*CW*(qup7*Mf2(4,Gen4)))*
     -              USfC(Sfe4,2,4,Gen4)
               qup27 = 3*CW*((qup13*MUE*TB + qup14*AfC(3,Gen4,Gen4))*
     -              Mf(3,Gen4)*USfC(Sfe5,2,3,Gen4)) + 
     -              qup17*TB*USfC(Sfe5,1,3,Gen4)
               qup28 = 3*CW*((qup15*MUEC*TB + qup16*Af(3,Gen4,Gen4))*
     -              Mf(3,Gen4)*USfC(Sfe5,1,3,Gen4)) + 
     -              2*TB*(2*MW*MZ*SW2*(qup9*SB) + 3*CW*(qup10*Mf2(3,Gen4)))*
     -              USfC(Sfe5,2,3,Gen4)
               qup29 = 3*CW*((qup12*MUE + qup15*AfC(4,Gen4,Gen4))*
     -              Mf(4,Gen4)*USfC(Sfe5,2,4,Gen4)) + 
     -              qup18*USfC(Sfe5,1,4,Gen4)
               qup30 = 3*CW*((qup11*MUEC + qup13*Af(4,Gen4,Gen4))*
     -              Mf(4,Gen4)*USfC(Sfe5,1,4,Gen4)) + 
     -              2*(MW*MZ*SW2*(CB*qup9) + 3*CW*(qup7*Mf2(4,Gen4)))*
     -              USfC(Sfe5,2,4,Gen4)

               loopa1 = loopa1 + Alfa/(CW2*MW*Pi*SW2)*(EL*(((1/36.D0*
     -              (qup21*Bval(bb0,iintQ6(Gen4,Sfe4,Sfe5))*(qup23*
     -              USf(Sfe5,1,3,Gen4) + qup24*USf(Sfe5,2,3,Gen4))) - 
     -              1/18.D0*(Cval(cc00,iintQ8(Gen4,Sfe4,Sfe5))*(qup19*
     -              (qup27*USf(Sfe4,1,3,Gen4) + qup28*
     -              USf(Sfe4,2,3,Gen4)) + qup21*(qup23*
     -              USf(Sfe5,1,3,Gen4) + qup24*USf(Sfe5,2,3,Gen4))))))/
     -              (SB*TB) + (-(1/72.D0*(qup22*
     -              Bval(bb0,iintQ5(Gen4,Sfe4,Sfe5))*(qup25*
     -              USf(Sfe5,1,4,Gen4) + qup26*USf(Sfe5,2,4,Gen4)))) + 
     -              (1/36.D0*(Cval(cc00,iintQ7(Gen4,Sfe4,Sfe5))))*
     -              (qup20*(qup29*USf(Sfe4,1,4,Gen4) + qup30*
     -              USf(Sfe4,2,4,Gen4)) + qup22*(qup25*
     -              USf(Sfe5,1,4,Gen4) +qup26*USf(Sfe5,2,4,Gen4))))/CB))

               loopa2 = loopa2 + Alfa/(CW2*MW*Pi*SW2)*(EL*((1/18.D0*
     -              (Cval(cc12,iintQ8(Gen4,Sfe4,Sfe5))*(qup19*(qup27*
     -              USf(Sfe4,1,3,Gen4) + qup28*USf(Sfe4,2,3,Gen4)) + 
     -              qup21*(qup23*USf(Sfe5,1,3,Gen4) + qup24*
     -              USf(Sfe5,2,3,Gen4)))))/(SB*TB) + ((-1/36.D0*
     -              (Cval(cc12,iintQ7(Gen4,Sfe4,Sfe5))))*(qup20*(qup29*
     -              USf(Sfe4,1,4,Gen4) + qup30*USf(Sfe4,2,4,Gen4)) + 
     -              qup22*(qup25*USf(Sfe5,1,4,Gen4) + 
     -              qup26*USf(Sfe5,2,4,Gen4))))/CB))
            enddo
         enddo
      enddo
      

** The lepton and slepton contributions
         if (sector .ne. 1) then
            lup1 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
            lup2 = -(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -           ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
            lup3 = cI*(SB*ZHggs(HiggsType,3))+SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)

            do Gen4 = 1,3
               loopa1 = loopa1 + Alfa/(CW*MW*Pi*SW2)*(EL*((1 - 4*SW2)/
     -              16.D0*(((q22+q12-q32)*lup2 + (q22+q12-q32)*lup3)*
     -              Cval(cc0,iintL1(Gen4))) + lup1*((-1 + 4*SW2)/4.D0*
     -              (Bval(bb0,iintL2(Gen4))) + (1 - 4*SW2)*
     -              (Cval(cc00,iintL1(Gen4)))))*Mf2(2,Gen4))/CB

               loopa2 = loopa2 + Alfa/(CW*MW*Pi*SW2)*(EL*((1 - 4*SW2)/
     -              16.D0*((2*lup2 + 2*lup3)*Cval(cc0,iintL1(Gen4))) + 
     -              lup1*((1 - 4*SW2)*
     -              (-Cval(cc12,iintL1(Gen4)))))*Mf2(2,Gen4))/CB

               loopa3 = loopa3 + Alfa/(CW*MW*Pi*SW2)*(EL*((1 - 4*SW2)/
     -              16.D0*((2*lup2 - 2*lup3)*Cval(cc0,iintL1(Gen4))))*
     -              Mf2(2,Gen4))/CB
            enddo

            lup4 = SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
            lup5 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
            lup6 = -(SAB*ZHggs(HiggsType,1)) + CAB*ZHggs(HiggsType,2)
            lup7 = -(cI*(CB*ZHggs(HiggsType,3))) + CA*
     -           ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2)
            lup8 = cI*(CB*ZHggs(HiggsType,3)) + CA*ZHggs(HiggsType,1) + 
     -           SA*ZHggs(HiggsType,2)
            lup9 = -(cI*(SB*ZHggs(HiggsType,3))) + SA*
     -           ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)
            lup10 = cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -           CA*ZHggs(HiggsType,2)

            do Sfe5 = 1,2
               do Sfe4 = 1,2
                  do Gen4 = 1,3
                     lup11 = MW*MZ*(-1 + 2*SW2)*(CB*lup5) + 2*CW*(lup4*
     -                    Mf2(2,Gen4))
                     lup12 = 2*SW2*(USf(Sfe5,1,2,Gen4)*
     -                    USfC(Sfe4,1,2,Gen4) + 
     -                    USf(Sfe5,2,2,Gen4)*USfC(Sfe4,2,2,Gen4)) - 
     -                    USf(Sfe5,1,2,Gen4)*USfC(Sfe4,1,2,Gen4)
                     lup13 = 2*SW2*(USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4) + 
     -                    USf(Sfe4,2,2,Gen4)*USfC(Sfe5,2,2,Gen4)) - 
     -                    USf(Sfe4,1,2,Gen4)*USfC(Sfe5,1,2,Gen4)
                     lup14 = CW*((lup8*MUE + lup10*AfC(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe4,2,2,Gen4)) + 
     -                    lup11*USfC(Sfe4,1,2,Gen4)
                     lup15 = CW*((lup7*MUEC + lup9*Af(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 2*(MW*MZ*
     -                    SW2*(CB*lup6) + CW*(lup4*Mf2(2,Gen4)))*
     -                    USfC(Sfe4,2,2,Gen4)
                     lup16 = CW*((lup8*MUE + lup10*AfC(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe5,2,2,Gen4)) + 
     -                    lup11*USfC(Sfe5,1,2,Gen4)
                     lup17 = CW*((lup7*MUEC + lup9*Af(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe5,1,2,Gen4)) + 2*(MW*MZ*
     -                    SW2*(CB*lup6) + CW*(lup4*Mf2(2,Gen4)))*
     -                    USfC(Sfe5,2,2,Gen4)

                     loopa1 = loopa1 + Alfa/(CW2*MW*Pi*SW2)*(EL*
     -                    ((-(1/8.D0*(lup13*
     -                    Bval(bb0,iintL3(Gen4,Sfe4,Sfe5))*(lup14*
     -                    USf(Sfe5,1,2,Gen4) + lup15*
     -                    USf(Sfe5,2,2,Gen4)))) + 1/4.D0*
     -                    (Cval(cc00,iintL4(Gen4,Sfe4,Sfe5))*(lup12*
     -                    (lup16*USf(Sfe4,1,2,Gen4) + lup17*
     -                    USf(Sfe4,2,2,Gen4)) + lup13*(lup14*
     -                    USf(Sfe5,1,2,Gen4) + lup15*
     -                    USf(Sfe5,2,2,Gen4)))))))/CB

                     loopa2 = loopa2 + Alfa/(CW2*MW*Pi*SW2)*(EL*
     -                    (-1/4.D0*(Cval(cc12,iintL4(Gen4,Sfe4,Sfe5))*
     -                    (lup12*(lup16*USf(Sfe4,1,2,Gen4) + 
     -                    lup17*USf(Sfe4,2,2,Gen4)) + 
     -                    lup13*(lup14*USf(Sfe5,1,2,Gen4) + 
     -                    lup15*USf(Sfe5,2,2,Gen4))))))/CB
                  enddo
               enddo
            enddo
         end if ! end of (s)lepton contributions


** the fermion contn to the form factors 
         if (fm .eq. 1) then
            a2vtxAZHfm = -loopa2         
            a1vtxAZHfm = loopa1 + (loopa2*q1q2) + ctFM         
            a3vtxAZHfm = loopa3*(0.D0,-1.D0)
         else
            a2vtxAZH2fm = -loopa2         
            a1vtxAZH2fm = loopa1 + (loopa2*q1q2) + ctFM         
            a3vtxAZH2fm = loopa3*(0.D0,-1.D0)
         end if

            
** And the SM bosons and their superpartners ...
         if (sector .ge. 3) then

            bup1 = SAB*ZHggs(HiggsType,1) - CAB*ZHggs(HiggsType,2)
            bup2 = SBA*ZHggs(HiggsType,1) + CBA*ZHggs(HiggsType,2)
            bup3 = 2*CW2*(SBA*ZHggs(HiggsType,1) + CBA*
     -           ZHggs(HiggsType,2)) + C2B*(SAB*ZHggs(HiggsType,1) - 
     -           CAB*ZHggs(HiggsType,2))

            loopa1 = loopa1+ (Alfa*MW)/Pi*(EL*(1/SW2*((Bval(bb0,iintB6)*
     -           (-((-CW2 + SW2)/(8.D0*CW*CW2)*(C2B*bup1)) + (-CW/2.D0 + 
     -           (-CW2 + SW2)/(4.D0*CW))*bup2) + Cval(cc00,iintB5)*
     -           ((-CW2 + SW2)/(2.D0*CW*CW2)*(C2B*bup1) + ((9*CW)/2.D0 - 
     -           (-2*CW2 + 4*SW2)/(4.D0*CW))*bup2)))+Cval(cc0,iintB5)*
     -           (1/(4.D0*CW)*(-((q22-q32)*bup2) - MW2/CW2*(C2B*bup1) + 
     -           (MW2*(-CW2 + SW2))/SW2*(bup2)) + 1/SW2*(bup2*((MW2*(2*
     -           CW2 + SW2))/(2.D0*CW) + CW*(5/4.D0*(q12+q22-q32) + 
     -           1/4.D0*(q12-q32))))) + 1/CW*((1/4.D0*(bup2*
     -           (Bval(bb0,iintB1) + Bval(bb0,iintB2))) + (-CW2 + SW2)/
     -           (CW2*SW2)*(bup3*(1/8.D0*Bval(bb0,iintB3) - 1/2.D0*
     -           Cval(cc00,iintB4)))) + 1/SW2*(1/4.D0*(bup2*
     -           (Cval(cc1,iintB5)*(SW2*(-(q12-q22-q32)) + CW2*((-q32 + 
     -           q12 + 3*q22) - (q12-q22-q32))) + Cval(cc2,iintB5)*
     -           (-SW2*(q12-q22+q32) + CW2*(3*(q12 - 1/3.D0*q32 + 
     -           1/3.D0*q22) - (q12-q22+q32)))))))))

            loopa2 = loopa2 +(Alfa*MW)/Pi*(EL*(1/SW2*(Cval(cc12,iintB5)*
     -           (-((-CW2 + SW2)/(2.D0*CW*CW2)*(C2B*bup1)) + 
     -           ((-9*CW)/2.D0 + (-2*CW2 + 4*SW2)/(4.D0*CW))*bup2))+
     -           Cval(cc0,iintB5)*(1/(4.D0*CW)*(-(4*bup2)) +1/SW2*(bup2*
     -           (CW*(5/4.D0*(8/5D0) + 1/4.D0*4)))) + 1/CW*(1/SW2*
     -           ((-CW2 + SW2)/(2.D0*CW2)*(bup3*Cval(cc12,iintB4)) + 
     -           1/4.D0*(bup2*(Cval(cc2,iintB5)*((-2*CW2 + 2*SW2) - 
     -           SW2*2 + CW2*(3*(-2/3D0) - (-4)))))))))

** the neutralino / chargino contribution
            if (sector .ge. 4) then
               bup8 = -(cI*(CB*ZHggs(HiggsType,3))) + CA*ZHggs(HiggsType,1) + 
     -              SA*ZHggs(HiggsType,2)
               bup9 = -(cI*(CA*ZHggs(HiggsType,1) + SA*ZHggs(HiggsType,2))) + 
     -              CB*ZHggs(HiggsType,3)
               bup10 = cI*(SB*ZHggs(HiggsType,3)) + SA*ZHggs(HiggsType,1) - 
     -              CA*ZHggs(HiggsType,2)
               bup11 = cI*(SA*ZHggs(HiggsType,1) - CA*ZHggs(HiggsType,2)) + 
     -              SB*ZHggs(HiggsType,3)

               do Cha5 = 1,2
                  do Cha4 = 1,2
                     bup4 = 2*SW2*IndexDelta(Cha4,Cha5) -2*UCha(Cha5,1)*
     -                    UChaC(Cha4,1) - UCha(Cha5,2)*UChaC(Cha4,2)
                     bup5 = 2*SW2*IndexDelta(Cha4,Cha5) -2*UCha(Cha4,1)*
     -                    UChaC(Cha5,1) - UCha(Cha4,2)*UChaC(Cha5,2)
                     bup6 = 2*SW2*IndexDelta(Cha4,Cha5) -2*VCha(Cha5,1)*
     -                    VChaC(Cha4,1) - VCha(Cha5,2)*VChaC(Cha4,2)
                     bup7 = 2*SW2*IndexDelta(Cha4,Cha5) -2*VCha(Cha4,1)*
     -                    VChaC(Cha5,1) - VCha(Cha4,2)*VChaC(Cha5,2)
                     bup12 =cI*(bup7*(bup10*UChaC(Cha5,2)*VChaC(Cha4,1)- 
     -                    bup8*UChaC(Cha5,1)*VChaC(Cha4,2)) + 
     -                    bup5*(bup10*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    bup8*UChaC(Cha4,1)*VChaC(Cha5,2))) + 
     -                    bup6*(bup11*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup9*UCha(Cha5,1)*VCha(Cha4,2)) + 
     -                    bup4*(bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2))
                     bup13 = cI*(bup4*
     -                    (bup10*UChaC(Cha5,2)*VChaC(Cha4,1) - 
     -                    bup8*UChaC(Cha5,1)*VChaC(Cha4,2)) + 
     -                    bup6*(bup10*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    bup8*UChaC(Cha4,1)*VChaC(Cha5,2))) + 
     -                    bup5*(bup11*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup9*UCha(Cha5,1)*VCha(Cha4,2)) + 
     -                    bup7*(bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2))

                     bup14 = cI*((q12+q22-q32)*
     -                    (bup4*(bup10*UChaC(Cha5,2)*VChaC(Cha4,1) - 
     -                    bup8*UChaC(Cha5,1)*VChaC(Cha4,2)) + 
     -                    bup6*(bup10*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    bup8*UChaC(Cha4,1)*VChaC(Cha5,2)))) + 
     -                    (q12+q22-q32)*(bup5*(bup11*UCha(Cha5,2)*
     -                    VCha(Cha4,1) + 
     -                    bup9*UCha(Cha5,1)*VCha(Cha4,2)) + 
     -                    bup7*(bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2)))

                     loopa1 = loopa1 + (Alfa*cI)/(CW*Pi*SW2)*(EL*
     -                    (1/8.D0*(bup14*Cval(cc0,iintB8(Cha4,Cha5))*
     -                    MCha(Cha4)) - 1/4.D0*(bup12*
     -                    Bval(bb0,iintB7(Cha4,Cha5))*MCha(Cha5)) + 
     -                    1/4.D0*(q12*Cval(cc1,iintB8(Cha4,Cha5))*
     -                    (bup13*MCha(Cha4) - bup12*MCha(Cha5))) + 
     -                    1/2.D0*((Cval(cc00,iintB8(Cha4,Cha5)))*
     -                    (bup13*MCha(Cha4) + bup12*MCha(Cha5))) + 
     -                    Cval(cc2,iintB8(Cha4,Cha5))*(1/8.D0*(bup14*
     -                    MCha(Cha4)) - 1/8.D0*(MCha(Cha5)*(cI*((q12+
     -                    q22-q32)*(bup7*(bup10*UChaC(Cha5,2)*
     -                    VChaC(Cha4,1) - bup8*UChaC(Cha5,1)*
     -                    VChaC(Cha4,2)) + bup5*(bup10*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup8*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2)))) + (q12+q22-q32)*(bup6*
     -                    (bup11*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup9*UCha(Cha5,1)*VCha(Cha4,2)) + bup4*
     -                    (bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2))))))))/
     -                    Sqrt(2.D0)

                     bup14 = cI*(2*(bup4*(bup10*UChaC(Cha5,2)*
     -                    VChaC(Cha4,1) - bup8*UChaC(Cha5,1)*
     -                    VChaC(Cha4,2)) + bup6*(bup10*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup8*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2)))) +2*(bup5*(bup11*UCha(Cha5,2)*
     -                    VCha(Cha4,1)+bup9*UCha(Cha5,1)*VCha(Cha4,2)) + 
     -                    bup7*(bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2)))
                     
                     loopa2 = loopa2 + (Alfa*cI)/(CW*Pi*SW2)*(EL*
     -                    (1/8.D0*(bup14*Cval(cc0,iintB8(Cha4,Cha5))*
     -                    MCha(Cha4)) + 
     -                    1/2.D0*((-Cval(cc12,iintB8(Cha4,Cha5)))*
     -                    (bup13*MCha(Cha4) + bup12*MCha(Cha5))) + 
     -                    Cval(cc2,iintB8(Cha4,Cha5))*
     -                    (1/8.D0*(bup14*MCha(Cha4))-1/8.D0*(MCha(Cha5)*
     -                    (cI*(2*(bup7*(bup10*UChaC(Cha5,2)*
     -                    VChaC(Cha4,1) - bup8*UChaC(Cha5,1)*
     -                    VChaC(Cha4,2)) + bup5*(bup10*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup8*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2)))) +2*(bup6*(bup11*UCha(Cha5,2)*
     -                    VCha(Cha4,1)+bup9*UCha(Cha5,1)*VCha(Cha4,2)) + 
     -                    bup4*(bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2))))))))/
     -                    Sqrt(2.D0)

                     bup14 = cI*(2*(bup4*(bup10*UChaC(Cha5,2)*
     -                    VChaC(Cha4,1) - 
     -                    bup8*UChaC(Cha5,1)*VChaC(Cha4,2)) + 
     -                    bup6*(bup10*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    bup8*UChaC(Cha4,1)*VChaC(Cha5,2)))) + 
     -                    (-2)*(bup5*(bup11*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    bup9*UCha(Cha5,1)*VCha(Cha4,2)) + 
     -                    bup7*(bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2)))

                     loopa3 = loopa3 + (Alfa*cI)/(CW*Pi*SW2)*(EL*
     -                    (1/8.D0*(bup14*Cval(cc0,iintB8(Cha4,Cha5))*
     -                    MCha(Cha4)) + Cval(cc2,iintB8(Cha4,Cha5))*
     -                    (1/8.D0*(bup14*MCha(Cha4))-1/8.D0*(MCha(Cha5)*
     -                    (cI*(2*(bup7*(bup10*UChaC(Cha5,2)*
     -                    VChaC(Cha4,1) - bup8*UChaC(Cha5,1)*
     -                    VChaC(Cha4,2)) + bup5*(bup10*UChaC(Cha4,2)*
     -                    VChaC(Cha5,1) - bup8*UChaC(Cha4,1)*
     -                    VChaC(Cha5,2)))) -2*(bup6*(bup11*UCha(Cha5,2)*
     -                    VCha(Cha4,1)+bup9*UCha(Cha5,1)*VCha(Cha4,2)) + 
     -                    bup4*(bup11*UCha(Cha4,2)*VCha(Cha5,1) + 
     -                    bup9*UCha(Cha4,1)*VCha(Cha5,2))))))))/
     -                    Sqrt(2.D0)
                  enddo
               enddo
            end if
         end if
** end of 'bosonic' contribution


** The eps part of the loop - the second part converts between FC notation
** and standard notation (in FC Eps = -i*levi-civita)
         loopa3 = loopa3*(0.D0,-1.D0)

         if (debug .ne. 0) then
            write(*,*)'photon-z-h loopa1 =', loopa1
            write(*,*)'photon-z-h loopa2 =', loopa2
            write(*,*)'photon-z-h loopa3 =', loopa3
         end if


** Now we want the form factors themselves
         a2vtxAZ = -loopa2         
         a1vtxAZ = loopa1 + (loopa2*q1q2) + ct          
         a3vtxAZ = loopa3
         
         if (debug .ne. 0) then
            write(*,*)"photon-z-higgs vtx a1 =", a1vtxAZ
            write(*,*)"photon-z-higgs vtx a2 =", a2vtxAZ
            write(*,*)"photon-z-higgs vtx a3 =", a3vtxAZ
            if (fm .eq. 1) then
               write(*,*)'FERMION contn to a1vtxAZ = ', a1vtxAZHfm
               write(*,*)'FERMION contn to a2vtxAZ = ', a2vtxAZHfm
               write(*,*)'FERMION contn to a3vtxAZ = ', a3vtxAZHfm
            else
               write(*,*)'FERMION contn to a1vtxAZ = ', a1vtxAZH2fm
               write(*,*)'FERMION contn to a2vtxAZ = ', a2vtxAZH2fm
               write(*,*)'FERMION contn to a3vtxAZ = ', a3vtxAZH2fm
            end if
         end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the WW self energy
*** in susy

      subroutine ww_selfSUSY(a1slfW,q12,q22,q32,q1q2,debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the z bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Integral arguments
      integer*8 iint1(3), iint2(3), iint3(3,2,2), iint4(3,2,2)
      double complex cint1(3), cint2(3,2), cint3(3,2)
      double complex cintL1(3), cintL2(3), cintL3(3,2)
      integer*8 iintL1(3), iintL2(3), iintL3(3,2), iintL4(3,2)
      double complex cintB1, cintB2, cintB3, cintB4, cintB5, cintB6
      double complex cintB7(4)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6, iintB7
      integer*8 iintB8, iintB9, iintB10, iintB11, iintB12, iintB13
      integer*8 iintB14, iintB15(2,4), iintB16(2,4)

** These are the contributions from the counter-term (g-mu-nu) and the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1


** This is the a1 formfactors for WW self energy.  
      double complex a1slfW


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


** Particle generation do loop parameters
      integer g4ini, Gen4, Sfe5, Sfe4, Cha4, Neu4

** Dummy variables
      double complex dup1
      double complex eta1, eta2, eta3, eta4, eta5


      double precision Alfa, Alfa2

      integer debug

** Parameter declarations and looptools
#include "looptools.h"

#include "mssm.inc"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions, 3 or 4 for SM-type
** particles and 5 or 6 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else 
         g4ini = 1
      end if


** Now for the arguments of the integrals.  Note that these are the same for
** heavy and light higgs.
** First the quarks and squarks
      do Gen4 = g4ini, 3
         cint1(Gen4) = A0(Mf2(4,Gen4))
         iint1(Gen4) = Bget(q22,Mf2(3,Gen4),Mf2(4,Gen4))
         iint2(Gen4) = Bget(q12,Mf2(3,Gen4),Mf2(4,Gen4))
      enddo

      do Sfe4 = 1,2
         do Gen4 = g4ini, 3
            cint2(Gen4,Sfe4) = A0(MSf2(Sfe4,3,Gen4))            
            cint3(Gen4,Sfe4) = A0(MSf2(Sfe4,4,Gen4))            
         enddo
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3               
               iint3(Gen4,Sfe4,Sfe5) = 
     -              Bget(q22,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,4,Gen4))
               iint4(Gen4,Sfe4,Sfe5) = 
     -              Bget(q12,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,4,Gen4))
            enddo
         enddo
      enddo

** Now the leptons and sleptons
      if (sector .ne. 1) then	
         do Gen4 = 1,3
            cintL1(Gen4) = A0(Mf2(2,Gen4))
            cintL2(Gen4) = A0(MSf2(1,1,Gen4))
            iintL1(Gen4) = Bget(q22,0.D0,Mf2(2,Gen4))
            iintL2(Gen4) = Bget(q12,0.D0,Mf2(2,Gen4))
         enddo

         do Sfe4 = 1,2
            do Gen4 = 1,3
               cintL3(Gen4,Sfe4) = A0(MSf2(Sfe4,2,Gen4))
               iintL3(Gen4,Sfe4) = Bget(q22,MSf2(1,1,Gen4),
     -              MSf2(Sfe4,2,Gen4))
               iintL4(Gen4,Sfe4) = Bget(q12,MSf2(1,1,Gen4),
     -              MSf2(Sfe4,2,Gen4))
            enddo
         enddo

      end if


** Integrals for SM bosons and their superpartners
      if (sector .ge. 3) then
         cintB1 = A0(MA02)
         cintB2 = A0(Mh02)
         cintB3 = A0(MHH2)
         cintB4 = A0(MHp2)
         cintB5 = A0(MW2)
         cintB6 = A0(MZ2)
         iintB1 = Bget(q22,0.D0,MW2)
         iintB2 = Bget(q22,Mh02,MW2)
         iintB3 = Bget(q22,MHH2,MW2)
         iintB4 = Bget(q22,MW2,MZ2)
         iintB5 = Bget(q22,MA02,MHp2)
         iintB6 = Bget(q22,MHH2,MHp2)
         iintB7 = Bget(q22,Mh02,MHp2)
         iintB8 = Bget(q12,0.D0,MW2)
         iintB9 = Bget(q12,Mh02,MW2)
         iintB10 = Bget(q12,MHH2,MW2)
         iintB11 = Bget(q12,MW2,MZ2)
         iintB12 = Bget(q12,MA02,MHp2)
         iintB13 = Bget(q12,MHH2,MHp2)
         iintB14 = Bget(q12,Mh02,MHp2)
         
** chargino / neutralino integrals
         if (sector .ge. 5) then
            do Neu4 = 1,4
               cintB7(Neu4) = A0(MNeu2(Neu4))
            enddo
            do Neu4 = 1,4
               do Cha4 = 1,2
                  iintB15(Cha4,Neu4) = Bget(q22,MCha2(Cha4),MNeu2(Neu4))
                  iintB16(Cha4,Neu4) = Bget(q12,MCha2(Cha4),MNeu2(Neu4))
               enddo
            enddo
         end if
      end if


** Clearing the loop values
      loopa1 = (0D0,0D0)  
      ct = (0D0,0D0) 


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
      ct = (-(MW*(q22*dZW1*EL*SBA)) + MW*(dMWsq1*EL*SBA) + 
     -     MW*MW2*(dZW1*EL*SBA))/((-MW2 + q22)*SW) + 
     -     (-(MW*(q12*dZW1*EL*SBA)) + MW*(dMWsq1*EL*SBA) + 
     -     MW*MW2*(dZW1*EL*SBA))/((-MW2 + q12)*SW)

** counterterm for the fermions
      ctFM = (-(MW*(q22*dZW1fm*EL*SBA)) + MW*(dMWsq1fm*EL*SBA) + 
     -     MW*MW2*(dZW1fm*EL*SBA))/((-MW2 + q22)*SW) + 
     -     (-(MW*(q12*dZW1fm*EL*SBA)) + MW*(dMWsq1fm*EL*SBA) + 
     -     MW*MW2*(dZW1fm*EL*SBA))/((-MW2 + q12)*SW)

      if (debug .ne. 0) then
         write(*,*)'WW ct=', ct
         write(*,*)'WW FERMION contn to ct=', ctFM
      end if


** Now the loop - first the quarks and squarks
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + 1/(4.D0*Pi*SW2)*(3*Alfa*MW*
     -        (q22*EL*SBA*Bval(bb1,iint1(Gen4))) + Alfa*MW*(EL*SBA*
     -        (-6*Bval(bb00,iint1(Gen4)) + 3*cint1(Gen4) + 
     -        3*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4))))/
     -        ((-MW2 + q22)*SW) + 1/(4.D0*Pi*SW2)*(3*Alfa*MW*
     -        (q12*EL*SBA*Bval(bb1,iint2(Gen4))) + Alfa*MW*(EL*SBA*
     -        (-6*Bval(bb00,iint2(Gen4)) + 3*cint1(Gen4) + 
     -        3*Bval(bb0,iint2(Gen4))*Mf2(3,Gen4))))/((-MW2 +q12)*SW)
      enddo

      do Sfe4 = 1,2
         do Gen4 = g4ini, 3
            dup1 = (Alfa*MW)/(8.D0*Pi*SW2)*
     -           (EL*SBA*(-3*cint2(Gen4,Sfe4)*USf(Sfe4,1,3,Gen4)*
     -           USfC(Sfe4,1,3,Gen4) - 3*cint3(Gen4,Sfe4)*
     -           USf(Sfe4,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)))/SW
            
            loopa1 = dup1/(-MW2 + q22) + dup1/(-MW2 + q12) + loopa1
         enddo
      enddo
      
      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3
               loopa1 = loopa1 + (3*Alfa*MW)/(2.D0*Pi*SW2)*
     -              (EL*SBA*Bval(bb00,iint3(Gen4,Sfe4,Sfe5))*
     -              USf(Sfe4,1,3,Gen4)*USf(Sfe5,1,4,Gen4)*
     -              USfC(Sfe4,1,3,Gen4)*USfC(Sfe5,1,4,Gen4))/
     -              ((-MW2 + q22)*SW) + (3*Alfa*MW)/(2.D0*Pi*SW2)*
     -              (EL*SBA*Bval(bb00,iint4(Gen4,Sfe4,Sfe5))*
     -              USf(Sfe4,1,3,Gen4)*USf(Sfe5,1,4,Gen4)*
     -              USfC(Sfe4,1,3,Gen4)*USfC(Sfe5,1,4,Gen4))/
     -              ((-MW2 + q12)*SW)
            enddo
         enddo
      enddo
      

** Now the leptons and sleptons
      if (sector .ne. 1) then
         do Gen4 = 1,3
            loopa1 = loopa1 + 1/(8.D0*Pi*SW2)*(2*Alfa*MW*
     -           (q22*EL*SBA*Bval(bb1,iintL1(Gen4))) + Alfa*MW*(EL*
     -           SBA*(-4*Bval(bb00,iintL1(Gen4)) + 2*cintL1(Gen4) - 
     -           cintL2(Gen4))))/((-MW2 + q22)*SW) + 
     -           1/(8.D0*Pi*SW2)*(2*Alfa*MW*(q12*
     -           EL*SBA*Bval(bb1,iintL2(Gen4))) + Alfa*MW*(EL*SBA*
     -           (-4*Bval(bb00,iintL2(Gen4)) + 2*cintL1(Gen4) - 
     -           cintL2(Gen4))))/((-MW2 + q12)*SW)
         enddo
         
         do Sfe4 = 1,2
            do Gen4 = 1,3
               loopa1 = loopa1 + (Alfa*MW)/(8.D0*Pi*SW2)*
     -              (EL*SBA*(4*Bval(bb00,iintL3(Gen4,Sfe4))*
     -              USf(Sfe4,1,2,Gen4)*USfC(Sfe4,1,2,Gen4) - 
     -              cintL3(Gen4,Sfe4)*USf(Sfe4,1,2,Gen4)*
     -              USfC(Sfe4,1,2,Gen4)))/((-MW2 + q22)*SW) + 
     -              (Alfa*MW)/(8.D0*Pi*SW2)*
     -              (EL*SBA*(4*Bval(bb00,iintL4(Gen4,Sfe4))*
     -              USf(Sfe4,1,2,Gen4)*USfC(Sfe4,1,2,Gen4) - 
     -              cintL3(Gen4,Sfe4)*USf(Sfe4,1,2,Gen4)*
     -              USfC(Sfe4,1,2,Gen4)))/((-MW2 + q12)*SW)
            enddo
         enddo         
      end if


** fermion form factor 
      a1selfWfm = loopa1 + ctFM


** 'Bosonic' contribution ...
      if (sector .ge. 3) then
         loopa1 = loopa1 + ((Alfa*MW)/(6.D0*Pi*SW2)*(q22*EL*SBA)/SW + 
     -        (Alfa*MW*MW2*(1 - CW2 - SW2))/(2.D0*Pi*SW2)*
     -        (EL*SBA)/SW + Bval(bb0,iintB4)*
     -        ((5*Alfa*CW2*MW)/(4.D0*Pi*SW2)*(q22*EL*SBA)/SW + 
     -        (Alfa*MW*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        (EL*SBA)/SW) - (Alfa*MW)/(16.D0*Pi*SW2)*
     -        ((cintB1 + cintB2 + cintB3)*EL*SBA)/SW - 
     -        (Alfa*MW)/(8.D0*Pi*SW2)*(cintB4*EL*SBA)/SW + 
     -        (Alfa*MW*(-7 + 4*SW2))/(8.D0*Pi*SW2)*(cintB5*EL*SBA)/
     -        SW - (Alfa*(1/4.D0 + CW2)*MW)/(4.D0*Pi*SW2)*
     -        (cintB6*EL*SBA)/SW - (Alfa*MW*MW2)/(4.D0*Pi*SW2)*
     -        (EL*SBA**3*Bval(bb0,iintB2))/SW - (Alfa*MW*MW2)/
     -        (4.D0*Pi*SW2)*(CBA**2*EL*SBA*Bval(bb0,iintB3))/SW + 
     -        (2*Alfa*MW)/Pi*(EL*SBA*Bval(bb00,iintB1))/SW + 
     -        (Alfa*(1/4.D0 + 2*CW2)*MW)/(Pi*SW2)*
     -        (EL*SBA*Bval(bb00,iintB4))/SW + (Alfa*MW)/(4.D0*Pi*SW2)*
     -        (EL*SBA*Bval(bb00,iintB5))/SW + (Alfa*MW)/(4.D0*Pi*SW2)*
     -        (EL*SBA**3*(Bval(bb00,iintB2) + Bval(bb00,iintB6)))/SW + 
     -        (Alfa*MW)/(4.D0*Pi*SW2)*(CBA**2*EL*SBA*
     -        (Bval(bb00,iintB3) + Bval(bb00,iintB7)))/SW + 
     -        (Alfa*MW)/(2.D0*Pi)*(q22*EL*SBA*Bval(bb1,iintB1))/SW + 
     -        (Alfa*CW2*MW)/(2.D0*Pi*SW2)*
     -        (q22*EL*SBA*Bval(bb1,iintB4))/SW + 
     -        1/(4.D0*Pi)*(Bval(bb0,iintB1)*(5*Alfa*MW*(q22*EL*SBA) - 
     -        Alfa*MW*MW2*(EL*SBA)))/SW)/(-MW2 + q22)
         
         loopa1 = loopa1 + ((Alfa*MW)/(6.D0*Pi*SW2)*(q12*EL*SBA)/SW + 
     -        (Alfa*MW*MW2*(1 - CW2 - SW2))/(2.D0*Pi*SW2)*
     -        (EL*SBA)/SW + Bval(bb0,iintB11)*
     -        ((5*Alfa*CW2*MW)/(4.D0*Pi*SW2)*(q12*EL*SBA)/SW + 
     -        (Alfa*MW*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        (EL*SBA)/SW) - (Alfa*MW)/(16.D0*Pi*SW2)*
     -        ((cintB1 + cintB2 + cintB3)*EL*SBA)/SW - 
     -        (Alfa*MW)/(8.D0*Pi*SW2)*(cintB4*EL*SBA)/SW + 
     -        (Alfa*MW*(-7 + 4*SW2))/(8.D0*Pi*SW2)*
     -        (cintB5*EL*SBA)/SW - (Alfa*(1/4.D0 + CW2)*MW)/
     -        (4.D0*Pi*SW2)*(cintB6*EL*SBA)/SW - (Alfa*MW*MW2)/
     -        (4.D0*Pi*SW2)*(CBA**2*EL*SBA*Bval(bb0,iintB10))/SW - 
     -        (Alfa*MW*MW2)/(4.D0*Pi*SW2)*(EL*SBA**3*Bval(bb0,iintB9))/
     -        SW + (Alfa*(1/4.D0 + 2*CW2)*MW)/(Pi*SW2)*
     -        (EL*SBA*Bval(bb00,iintB11))/SW + (Alfa*MW)/(4.D0*Pi*SW2)*
     -        (EL*SBA*Bval(bb00,iintB12))/SW + 
     -        (Alfa*MW)/(4.D0*Pi*SW2)*(CBA**2*EL*SBA*
     -        (Bval(bb00,iintB10) + Bval(bb00,iintB14)))/SW + 
     -        (2*Alfa*MW)/Pi*(EL*SBA*Bval(bb00,iintB8))/SW + 
     -        (Alfa*MW)/(4.D0*Pi*SW2)*(EL*SBA**3*
     -        (Bval(bb00,iintB13) + Bval(bb00,iintB9)))/SW + 
     -        (Alfa*CW2*MW)/(2.D0*Pi*SW2)*(q12*EL*SBA*
     -        Bval(bb1,iintB11))/SW + (Alfa*MW)/(2.D0*Pi)*(q12*EL*
     -        SBA*Bval(bb1,iintB8))/SW + 1/(4.D0*Pi)*(Bval(bb0,iintB8)*
     -        (5*Alfa*MW*(q12*EL*SBA) - 
     -        Alfa*MW*MW2*(EL*SBA)))/SW)/(-MW2 + q12)


* charginos and neutralinos
         if (sector .ge. 5) then
            do Neu4 = 1,4
               do Cha4 = 1,2     
                  eta1 = (Sqrt(2.D0)*UChaC(Cha4,1)*ZNeu(Neu4,2) + 
     -                 UChaC(Cha4,2)*ZNeu(Neu4,3))/Sqrt(2.D0)
                  eta2 = (Sqrt(2.D0)*VChaC(Cha4,1)*ZNeu(Neu4,2) - 
     -                 VChaC(Cha4,2)*ZNeu(Neu4,4))/Sqrt(2.D0)
                  eta3 = (Sqrt(2.D0)*UCha(Cha4,1)*ZNeuC(Neu4,2) + 
     -                 UCha(Cha4,2)*ZNeuC(Neu4,3))/Sqrt(2.D0)
                  eta4 = (Sqrt(2.D0)*VCha(Cha4,1)*ZNeuC(Neu4,2) - 
     -                 VCha(Cha4,2)*ZNeuC(Neu4,4))/Sqrt(2.D0)
                  eta5 = Alfa/(2.D0*Pi*SW2)*
     -                 (MW*(eta2*eta4*EL*SBA*MCha2(Cha4)) + 
     -                 eta1*(MW*(eta3*EL*SBA*MCha2(Cha4)) - 
     -                 MW*(eta2*EL*SBA*MCha(Cha4)*MNeu(Neu4))) - 
     -                 MW*(eta3*eta4*EL*SBA*MCha(Cha4)*MNeu(Neu4)))/SW
                  
                  loopa1 = loopa1 + 1/(2.D0*Pi*SW2)*(2*Pi*SW2*
     -                 (eta5*SW*Bval(bb0,iintB15(Cha4,Neu4))) + 
     -                 Alfa*(q22*(MW*(eta1*eta3*EL*SBA*
     -                 Bval(bb1,iintB15(Cha4,Neu4)))+MW*(eta2*eta4*EL*SBA*
     -                 Bval(bb1,iintB15(Cha4,Neu4))))) + 
     -                 Alfa*(MW*(eta1*eta3*EL*SBA*(-2*
     -                 Bval(bb00,iintB15(Cha4,Neu4)) + cintB7(Neu4))) +
     -                 MW*(eta2*eta4*EL*SBA*
     -                 (-2*Bval(bb00,iintB15(Cha4,Neu4)) + cintB7(Neu4)))
     -                 ))/((-MW2 + q22)*SW) + 1/(2.D0*Pi*SW2)*(2*Pi*SW2*
     -                 (eta5*SW*Bval(bb0,iintB16(Cha4,Neu4))) + Alfa*
     -                 (q12*(MW*(eta1*eta3*EL*SBA*
     -                 Bval(bb1,iintB16(Cha4,Neu4)))+MW*(eta2*eta4*EL*SBA*
     -                 Bval(bb1,iintB16(Cha4,Neu4)))))+Alfa*(MW*(eta1*eta3*
     -                 EL*SBA*(-2*Bval(bb00,iintB16(Cha4,Neu4)) + 
     -                 cintB7(Neu4)))+ MW*(eta2*eta4*EL*SBA*
     -                 (-2*Bval(bb00,iintB16(Cha4,Neu4)) + cintB7(Neu4)))
     -                 ))/((-MW2 + q12)*SW)
               enddo
            enddo
         end if
      end if


      if (debug .ne. 0) then
         write(*,*)'WW loopa1=', loopa1
      end if


** Now we want the form factor itself         
      a1slfW = loopa1 + ct 
      
      if (debug .ne. 0) then
         write(*,*)"WW slf a1 =", a1slfW
         write(*,*)"WW FERMION contn to slf a1 =", a1selfWfm
      end if

      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the ZZ self energy
*** in susy

      subroutine zz_selfSUSY(a1slfZ,q12,q22,q32,q1q2,debug)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Dummy variables used when calculating the loop elements
      double complex dup1, dup2, dup3, dup4, dup5, dup6
      double complex eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8

** The arguments of the integrals
      double complex cint1(3), cint2(3), cint3(3,2), cint4(3,2)
      integer*8 iint1(3), iint2(3), iint3(3), iint4(3), iint5(3,2,2)
      integer*8 iint6(3,2,2), iint7(3,2,2), iint8(3,2,2)
      double complex cintL1(3), cintL2(3), cintL3(3,2)
      integer*8 iintL1, iintL2, iintL3(3), iintL4(3), iintL5(3), iintL6(3)
      integer*8 iintL7(3,2,2), iintL8(3,2,2)
      double complex cintB1, cintB2, cintB3, cintB4, cintB5, cintB6
      double complex cintB7(4), cintB8(2)
      integer*8 iintB1, iintB2, iintB3, iintB4, iintB5, iintB6, iintB7
      integer*8 iintB8, iintB9, iintB10, iintB11, iintB12, iintB13(4,4)
      integer*8 iintB14(4,4), iintB15(2,2), iintB16(2,2)

** These are the contributions from the counter-term (g-mu-nu), the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1


** This is the formfactor (a1: g-mu-nu) for ZZ-higgs.  
      double complex a1slfZ


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


** Particle generation do loop variables
      integer Gen4, Sfe4, Sfe5, g4ini
      integer Neu4, Neu5, Cha4, Cha5

      double precision Alfa, Alfa2
      double precision IndexDelta(4,4)

      integer debug

** Model parameters
#include "mssm.inc"   
#include "global.inc"   

** This is needed for looptools. 
#include "looptools.h"


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions and 3 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else 
         g4ini = 1
      end if


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      do Neu4 = 1, 4
         do Neu5 = 1, 4
            if (Neu4 .eq. Neu5) then
               IndexDelta(Neu4,Neu5) = 1D0
            else
               IndexDelta(Neu4,Neu5) = 0D0
            end if
         end do
      end do


** Now for the arguments of the integrals.  Note that these are the same
** for both the heavy and light higgs contributions.
** First the quarks and squarks
      do Gen4 = g4ini, 3         
         cint1(Gen4) = A0(Mf2(3,Gen4))         
         cint2(Gen4) = A0(Mf2(4,Gen4))         
         iint1(Gen4) = Bget(q22,Mf2(3,Gen4),Mf2(3,Gen4))         
         iint2(Gen4) = Bget(q22,Mf2(4,Gen4),Mf2(4,Gen4))         
         iint3(Gen4) = Bget(q12,Mf2(3,Gen4),Mf2(3,Gen4))         
         iint4(Gen4) = Bget(q12,Mf2(4,Gen4),Mf2(4,Gen4))         
      enddo

      do Sfe4 = 1,2
         do Gen4 = g4ini, 3
            cint3(Gen4,Sfe4) = A0(MSf2(Sfe4,3,Gen4))
            cint4(Gen4,Sfe4) = A0(MSf2(Sfe4,4,Gen4))
         enddo
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3
               iint5(Gen4,Sfe4,Sfe5) = 
     -              Bget(q22,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,3,Gen4))
               iint6(Gen4,Sfe4,Sfe5) = 
     -              Bget(q22,MSf2(Sfe4,4,Gen4),MSf2(Sfe5,4,Gen4))
               iint7(Gen4,Sfe4,Sfe5) = 
     -              Bget(q12,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,3,Gen4))
               iint8(Gen4,Sfe4,Sfe5) = 
     -              Bget(q12,MSf2(Sfe4,4,Gen4),MSf2(Sfe5,4,Gen4))
            enddo
         enddo
      enddo

** Now the leptons and sleptons
      if (sector .ne. 1) then
         iintL1 = Bget(q22,0.D0,0.D0)
         iintL2 = Bget(q12,0.D0,0.D0)

         do Gen4 = 1,3
            cintL1(Gen4) = A0(Mf2(2,Gen4))
            cintL2(Gen4) = A0(MSf2(1,1,Gen4))
            iintL3(Gen4) = Bget(q22,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL4(Gen4) = Bget(q22,MSf2(1,1,Gen4),MSf2(1,1,Gen4))
            iintL5(Gen4) = Bget(q12,Mf2(2,Gen4),Mf2(2,Gen4))
            iintL6(Gen4) = Bget(q12,MSf2(1,1,Gen4),MSf2(1,1,Gen4))
         enddo

         do Sfe4 = 1,2
            do Gen4 = 1,3
               cintL3(Gen4,Sfe4) = A0(MSf2(Sfe4,2,Gen4))
            enddo
         enddo

         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  iintL7(Gen4,Sfe4,Sfe5) = 
     -                 Bget(q22,MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4))
                  iintL8(Gen4,Sfe4,Sfe5) = 
     -                 Bget(q12,MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4))
               enddo
            enddo
         enddo
      end if

** Integrals for SM bosons and their superpartners ...
      if (sector .ge. 3) then
         cintB1 = A0(MHp2)
         cintB2 = A0(MW2)
         cintB3 = A0(MA02)
         cintB4 = A0(Mh02)
         cintB5 = A0(MHH2)
         cintB6 = A0(MZ2)
         iintB1 = Bget(q22,Mh02,MZ2)
         iintB2 = Bget(q22,MHH2,MZ2)
         iintB3 = Bget(q22,MW2,MW2)
         iintB4 = Bget(q22,MA02,MHH2)
         iintB5 = Bget(q22,MA02,Mh02)
         iintB6 = Bget(q22,MHp2,MHp2)
         iintB7 = Bget(q12,Mh02,MZ2)
         iintB8 = Bget(q12,MHH2,MZ2)
         iintB9 = Bget(q12,MW2,MW2)
         iintB10 = Bget(q12,MA02,MHH2)
         iintB11 = Bget(q12,MA02,Mh02)
         iintB12 = Bget(q12,MHp2,MHp2)

** neutralino / chargino integrals
         if (sector .ge. 5) then
            do Cha5 = 1,2
               cintB8(Cha5) = A0(MCha2(Cha5))
            enddo
            do Neu5 = 1,4
               cintB7(Neu5) = A0(MNeu2(Neu5))
            enddo
            do Cha5 = 1,2
               do Cha4 = 1,2
                  iintB15(Cha4,Cha5) = Bget(q22,MCha2(Cha4),MCha2(Cha5))
                  iintB16(Cha4,Cha5) = Bget(q12,MCha2(Cha4),MCha2(Cha5))
               enddo
            enddo
            do Neu5 = 1,4
               do Neu4 = 1,4
                  iintB13(Neu4,Neu5) = Bget(q22,MNeu2(Neu4),MNeu2(Neu5))
                  iintB14(Neu4,Neu5) = Bget(q12,MNeu2(Neu4),MNeu2(Neu5))
               enddo
            enddo
         end if
         
      end if


** Clearing the values of each contribution
      loopa1 = (0D0,0D0)  
      ct = (0D0,0D0)   


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
      ct = 1/CW2*(-(MW*(q22*dZZZ1*EL*SBA)) +  MW*(dMZsq1*
     -     EL*SBA) + MW*MZ2*(dZZZ1*EL*SBA))/((-MZ2 + q22)*SW) + 
     -     1/CW2*(-(MW*(q12*dZZZ1*EL*SBA)) + MW*(dMZsq1*EL*
     -     SBA) + MW*MZ2*(dZZZ1*EL*SBA))/((-MZ2 + q12)*SW)

** fermion value of the counterterm
      ctFM = 1/CW2*(-(MW*(q22*dZZZ1fm*EL*SBA)) +  MW*(dMZsq1fm*
     -     EL*SBA) + MW*MZ2*(dZZZ1fm*EL*SBA))/((-MZ2 + q22)*SW) + 
     -     1/CW2*(-(MW*(q12*dZZZ1fm*EL*SBA)) + MW*(dMZsq1fm*EL*
     -     SBA) + MW*MZ2*(dZZZ1fm*EL*SBA))/((-MZ2 + q12)*SW)

      if (debug .ne. 0) then
         write(*,*)'ZZ ct=', ct
         write(*,*)'ZZ FERMION contn to the ct=', ctFM
      end if
         

** Now the loop: first the quarks and squarks
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + ((Alfa*MW*(-9 + (24 - 32*SW2)*SW2))/
     -        (12.D0*CW2**2*Pi*SW2)*(EL*SBA*Bval(bb00,iint1(Gen4)))/
     -        SW + (Alfa*MW*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2**2*
     -        Pi*SW2)*(EL*SBA*Bval(bb00,iint2(Gen4)))/SW + (Alfa*
     -        MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q22*EL*SBA*Bval(bb1,iint1(Gen4)))/SW + (Alfa*
     -        MW*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q22*EL*SBA*Bval(bb1,iint2(Gen4)))/SW + (Alfa*
     -        MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (EL*SBA*cint1(Gen4))/SW + (Alfa*MW*(9 + 
     -        SW2*(-12 + 8*SW2)))/(24.D0*CW2**2*Pi*SW2)*(EL*SBA*
     -        cint2(Gen4))/SW + (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -        (EL*SBA*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(EL*SBA*
     -        Bval(bb0,iint2(Gen4))*Mf2(4,Gen4))/SW)/(-MZ2 + q22)
         
         loopa1 = loopa1 + ((Alfa*MW*(-9 + (24 - 32*SW2)*SW2))/
     -        (12.D0*CW2**2*Pi*SW2)*(EL*SBA*Bval(bb00,iint3(Gen4)))/
     -        SW + (Alfa*MW*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2**2*
     -        Pi*SW2)*(EL*SBA*Bval(bb00,iint4(Gen4)))/SW + (Alfa*
     -        MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q12*EL*SBA*Bval(bb1,iint3(Gen4)))/SW + (Alfa*
     -        MW*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (q12*EL*SBA*Bval(bb1,iint4(Gen4)))/SW + (Alfa*
     -        MW*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2**2*Pi*SW2)*
     -        (EL*SBA*cint1(Gen4))/SW + (Alfa*MW*(9 + SW2*(-12 + 
     -        8*SW2)))/(24.D0*CW2**2*Pi*SW2)*(EL*SBA*cint2(Gen4))/SW+ 
     -        (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -        (EL*SBA*Bval(bb0,iint3(Gen4))*Mf2(3,Gen4))/SW + 
     -        (3*Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(EL*SBA*
     -        Bval(bb0,iint4(Gen4))*Mf2(4,Gen4))/SW)/(-MZ2 + q12)
      enddo
         
      do Sfe4 = 1,2
         do Gen4 = g4ini, 3
            dup1 = (3 - 4*SW2)**2*
     -           (USf(Sfe4,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)) + 
     -           16*SW2**2*(USf(Sfe4,2,3,Gen4)*USfC(Sfe4,2,3,Gen4))
            dup2 = (3 - 2*SW2)**2*
     -           (USf(Sfe4,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)) + 
     -           4*SW2**2*(USf(Sfe4,2,4,Gen4)*USfC(Sfe4,2,4,Gen4))
            
            loopa1 = loopa1 + Alfa/(24.D0*CW2**2*Pi*SW2)*
     -           (-(MW*(dup1*EL*SBA*cint3(Gen4,Sfe4))) - 
     -           MW*(dup2*EL*SBA*cint4(Gen4,Sfe4)))/((-MZ2+q22)*SW) + 
     -           Alfa/(24.D0*CW2**2*Pi*SW2)*
     -           (-(MW*(dup1*EL*SBA*cint3(Gen4,Sfe4))) - 
     -           MW*(dup2*EL*SBA*cint4(Gen4,Sfe4)))/((-MZ2 + q12)*SW)
         enddo
      enddo
      
      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3
               dup3 = SW2*(4*USf(Sfe5,1,3,Gen4)*USfC(Sfe4,1,3,Gen4) + 
     -              4*USf(Sfe5,2,3,Gen4)*USfC(Sfe4,2,3,Gen4)) - 
     -              3*USf(Sfe5,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)
               dup4 = SW2*(2*USf(Sfe5,1,4,Gen4)*USfC(Sfe4,1,4,Gen4) + 
     -              2*USf(Sfe5,2,4,Gen4)*USfC(Sfe4,2,4,Gen4)) - 
     -              3*USf(Sfe5,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)
               dup5 = SW2*(4*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4) + 
     -              4*USf(Sfe4,2,3,Gen4)*USfC(Sfe5,2,3,Gen4)) - 
     -              3*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)
               dup6 = SW2*(2*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4) + 
     -              2*USf(Sfe4,2,4,Gen4)*USfC(Sfe5,2,4,Gen4)) - 
     -              3*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)
               
               loopa1 = loopa1+ Alfa/(12.D0*CW2**2*Pi*SW2)*(MW*(dup3*
     -              dup5*EL*SBA*Bval(bb00,iint5(Gen4,Sfe4,Sfe5))) + 
     -              MW*(dup4*dup6*EL*SBA*
     -              Bval(bb00,iint6(Gen4,Sfe4,Sfe5))))/
     -              ((-MZ2 + q22)*SW) + Alfa/(12.D0*CW2**2*Pi*SW2)*
     -              (MW*(dup3*dup5*EL*SBA*
     -              Bval(bb00,iint7(Gen4,Sfe4,Sfe5))) + 
     -              MW*(dup4*dup6*EL*SBA*
     -              Bval(bb00,iint8(Gen4,Sfe4,Sfe5))))/((-MZ2+q12)*SW)
            enddo
         enddo
      enddo
      

** Now the leptons and sleptons
      if (sector .ne. 1) then
            loopa1 = loopa1 + 1/(8.D0*CW2**2*Pi*SW2)*
     -           (-(6*Alfa*MW*(EL*SBA*Bval(bb00,iintL1))) + 
     -           3*Alfa*MW*(q22*EL*SBA*Bval(bb1,iintL1)))/((-MZ2 +  
     -           q22)*SW) + 1/(8.D0*CW2**2*Pi*SW2)*
     -           (-(6*Alfa*MW*(EL*SBA*Bval(bb00,iintL2))) + 3*Alfa*MW*
     -           (q12*EL*SBA*Bval(bb1,iintL2)))/((-MZ2 + q12)*SW)

            do Gen4 = 1,3
               loopa1 = loopa1 + ((Alfa*MW*(-1 + (4 - 8*SW2)*SW2))/
     -              (4.D0*CW2**2*Pi*SW2)*
     -              (EL*SBA*Bval(bb00,iintL3(Gen4)))/SW + 
     -              (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -              (EL*SBA*Bval(bb00,iintL4(Gen4)))/SW + (Alfa*MW*
     -              (1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -              (q22*EL*SBA*Bval(bb1,iintL3(Gen4)))/SW + (Alfa*MW*
     -              (1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -              (EL*SBA*cintL1(Gen4))/SW - (Alfa*MW)/(8.D0*CW2**
     -              2*Pi*SW2)*(EL*SBA*cintL2(Gen4))/SW + 
     -              (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -              (EL*SBA*Bval(bb0,iintL3(Gen4))*Mf2(2,Gen4))/SW)/
     -              (-MZ2 + q22) + ((Alfa*MW*(-1 + 
     -              (4 - 8*SW2)*SW2))/(4.D0*CW2**2*Pi*SW2)*
     -              (EL*SBA*Bval(bb00,iintL5(Gen4)))/SW + 
     -              (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -              (EL*SBA*Bval(bb00,iintL6(Gen4)))/SW + (Alfa*MW*
     -              (1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -              (q12*EL*SBA*Bval(bb1,iintL5(Gen4)))/SW + (Alfa*MW*
     -              (1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2**2*Pi*SW2)*
     -              (EL*SBA*cintL1(Gen4))/SW - (Alfa*MW)/(8.D0*CW2**2*
     -              Pi*SW2)*(EL*SBA*cintL2(Gen4))/SW + 
     -              (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(EL*SBA*
     -              Bval(bb0,iintL5(Gen4))*Mf2(2,Gen4))/SW)/(-MZ2 + q12)
            enddo
            
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  dup1 = (1 - 2*SW2)**2*
     -                 (USf(Sfe4,1,2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 
     -                 4*SW2**2*(USf(Sfe4,2,2,Gen4)*USfC(Sfe4,2,2,Gen4))
                  
                  loopa1 = loopa1- (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(dup1*
     -                 EL*SBA*cintL3(Gen4,Sfe4))/((-MZ2 + q22)*SW)-
     -                 (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -                 (dup1*EL*SBA*cintL3(Gen4,Sfe4))/((-MZ2 + q12)*SW)
               enddo
            enddo
            
            do Sfe5 = 1,2
               do Sfe4 = 1,2
                  do Gen4 = 1,3
                     dup2 = SW2*(2*USf(Sfe5,1,2,Gen4)*
     -                    USfC(Sfe4,1,2,Gen4) + 
     -                    2*USf(Sfe5,2,2,Gen4)*USfC(Sfe4,2,2,Gen4)) - 
     -                    USf(Sfe5,1,2,Gen4)*USfC(Sfe4,1,2,Gen4)
                     dup3 = SW2*(2*USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4) + 
     -                    2*USf(Sfe4,2,2,Gen4)*USfC(Sfe5,2,2,Gen4)) - 
     -                    USf(Sfe4,1,2,Gen4)*USfC(Sfe5,1,2,Gen4)
                     
                     loopa1 = loopa1 + (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -                    (dup2*dup3*EL*SBA*
     -                    Bval(bb00,iintL7(Gen4,Sfe4,Sfe5)))/((-MZ2 + 
     -                    q22)*SW)+(Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -                    (dup2*dup3*EL*SBA*
     -                    Bval(bb00,iintL8(Gen4,Sfe4,Sfe5)))/ 
     -                    ((-MZ2 + q12)*SW)
                  enddo
               enddo
            enddo
         
      end if


** The fermion contn to the form factors 
      a1selfZfm = loopa1 + ctFM


** And the sm bosons and their superpartners ...
      if (sector .ge. 3) then
         loopa1 = loopa1 + ((Alfa*MW)/(6.D0*Pi*SW2)*(q22*EL*SBA)/SW + 
     -        Bval(bb0,iintB3)*
     -        ((5*Alfa*MW)/(4.D0*Pi*SW2)*(q22*EL*SBA)/SW + 
     -        (Alfa*MW*MW2*(CW2**2 - SW2**2))/(2.D0*CW2**2*Pi*SW2)*
     -        (EL*SBA)/SW) + (Alfa*MW*(-SW2**2 + CW2*(-CW2 + 2*SW2)))/
     -        (8.D0*CW2**2*Pi*SW2)*(cintB1*EL*SBA)/SW + 
     -        (Alfa*MW*(-SW2**2 + CW2*(-9*CW2 + 2*SW2)))/
     -        (8.D0*CW2**2*Pi*SW2)*(cintB2*EL*SBA)/SW - 
     -        (Alfa*MW)/(16.D0*CW2**2*Pi*SW2)*
     -        ((cintB3 + cintB4 + cintB5 + cintB6)*EL*SBA)/SW - 
     -        (Alfa*MW*MW2)/(4.D0*CW2**3*Pi*SW2)*(EL*SBA**3*
     -        Bval(bb0,iintB1))/SW - (Alfa*MW*MW2)/(4.D0*CW2**3*Pi*SW2)*
     -        (CBA**2*EL*SBA*Bval(bb0,iintB2))/SW + 
     -        (Alfa*MW*(CW2*(9*CW2 - 2*SW2) + SW2**2))/
     -        (4.D0*CW2**2*Pi*SW2)*(EL*SBA*Bval(bb00,iintB3))/SW +
     -        (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*(EL*SBA**3*
     -        (Bval(bb00,iintB1) + Bval(bb00,iintB4)))/SW + 
     -        (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*(CBA**2*EL*SBA*
     -        (Bval(bb00,iintB2) + Bval(bb00,iintB5)))/SW + 
     -        (Alfa*MW*(CW2*(CW2 - 2*SW2) + SW2**2))/
     -        (4.D0*CW2**2*Pi*SW2)*(EL*SBA*Bval(bb00,iintB6))/SW +
     -        (Alfa*MW)/(2.D0*Pi*SW2)*
     -        (q22*EL*SBA*Bval(bb1,iintB3))/SW)/(-MZ2 + q22)
         
         loopa1 = loopa1 + ((Alfa*MW)/(6.D0*Pi*SW2)*(q12*EL*SBA)/SW + 
     -        Bval(bb0,iintB9)*((5*Alfa*MW)/(4.D0*Pi*SW2)*(q12*EL*SBA)/
     -        SW + (Alfa*MW*MW2*(CW2**2 - SW2**2))/(2.D0*CW2**2*Pi*SW2)*
     -        (EL*SBA)/SW) + (Alfa*MW*(-SW2**2 + CW2*(-CW2 + 2*SW2)))/
     -        (8.D0*CW2**2*Pi*SW2)*(cintB1*EL*SBA)/SW + 
     -        (Alfa*MW*(-SW2**2 + CW2*(-9*CW2 + 2*SW2)))/
     -        (8.D0*CW2**2*Pi*SW2)*(cintB2*EL*SBA)/SW - 
     -        (Alfa*MW)/(16.D0*CW2**2*Pi*SW2)*
     -        ((cintB3 + cintB4 + cintB5 + cintB6)*EL*SBA)/SW - 
     -        (Alfa*MW*MW2)/(4.D0*CW2**3*Pi*SW2)*
     -        (EL*SBA**3*Bval(bb0,iintB7))/SW - 
     -        (Alfa*MW*MW2)/(4.D0*CW2**3*Pi*SW2)*
     -        (CBA**2*EL*SBA*Bval(bb0,iintB8))/SW + 
     -        (Alfa*MW*(CW2*(CW2 -2*SW2)+ SW2**2))/(4.D0*CW2**2*Pi*SW2)*
     -        (EL*SBA*Bval(bb00,iintB12))/SW + 
     -        (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*(EL*SBA**3*
     -        (Bval(bb00,iintB10) + Bval(bb00,iintB7)))/SW + 
     -        (Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*(CBA**2*EL*SBA*
     -        (Bval(bb00,iintB11) + Bval(bb00,iintB8)))/SW + 
     -        (Alfa*MW*(CW2*(9*CW2 - 2*SW2) + SW2**2))/
     -        (4.D0*CW2**2*Pi*SW2)*(EL*SBA*Bval(bb00,iintB9))/SW +
     -        (Alfa*MW)/(2.D0*Pi*SW2)*
     -        (q12*EL*SBA*Bval(bb1,iintB9))/SW)/(-MZ2 + q12)
         
         
* charginos and neutralinos
         if (sector .ge. 5) then
            do Cha5 = 1,2
               do Cha4 = 1,2
                  eta4 = SW2*IndexDelta(Cha4,Cha5) - 
     -                 1/2.D0*(UCha(Cha5,2)*UChaC(Cha4,2)) - 
     -                 UCha(Cha5,1)*UChaC(Cha4,1)
                  eta5 = SW2*IndexDelta(Cha4,Cha5) - 
     -                 1/2.D0*(UCha(Cha4,2)*UChaC(Cha5,2)) - 
     -                 UCha(Cha4,1)*UChaC(Cha5,1)
                  eta6 = SW2*IndexDelta(Cha4,Cha5) - 
     -                 1/2.D0*(VCha(Cha5,2)*VChaC(Cha4,2)) - 
     -                 VCha(Cha5,1)*VChaC(Cha4,1)
                  eta7 = SW2*IndexDelta(Cha4,Cha5) - 
     -                 1/2.D0*(VCha(Cha4,2)*VChaC(Cha5,2)) - 
     -                 VCha(Cha4,1)*VChaC(Cha5,1)
                  eta8 = Alfa/(2.D0*CW2**2*Pi*SW2)*
     -                 (-(MW*(eta5*eta7*EL*SBA*MCha(Cha4)*MCha(Cha5))) + 
     -                 eta4*(-(MW*(eta6*EL*SBA*MCha(Cha4)*
     -                 MCha(Cha5))) + MW*(eta5*EL*SBA*MCha2(Cha4))) + 
     -                 MW*(eta6*eta7*EL*SBA*MCha2(Cha4)))/SW
                  
                  loopa1 = loopa1 + 1/(2.D0*CW2**2*Pi*SW2)*(2*CW2**2*Pi*
     -                 SW2*(eta8*SW*Bval(bb0,iintB15(Cha4,Cha5))) + 
     -                 Alfa*(q22*(MW*(eta4*eta5*EL*SBA*
     -                 Bval(bb1,iintB15(Cha4,Cha5)))+MW*(eta6*eta7*EL*SBA*
     -                 Bval(bb1,iintB15(Cha4,Cha5))))) + 
     -                 Alfa*(MW*(eta4*eta5*EL*SBA*
     -                 (-2*Bval(bb00,iintB15(Cha4,Cha5)) + cintB8(Cha5)))
     -                 + MW*(eta6*eta7*EL*SBA*
     -                 (-2*Bval(bb00,iintB15(Cha4,Cha5)) + cintB8(Cha5)))
     -                 ))/((-MZ2 + q22)*SW) + 1/(2.D0*CW2**2*Pi*SW2)*
     -                 (2*CW2**2*Pi*SW2*(eta8*SW*
     -                 Bval(bb0,iintB16(Cha4,Cha5))) + Alfa*(q12*(MW*
     -                 (eta4*eta5*EL*SBA*Bval(bb1,iintB16(Cha4,Cha5)))+
     -                 MW*(eta6*eta7*EL*SBA*
     -                 Bval(bb1,iintB16(Cha4,Cha5))))) + 
     -                 Alfa*(MW*(eta4*eta5*EL*SBA*
     -                 (-2*Bval(bb00,iintB16(Cha4,Cha5)) + cintB8(Cha5)))
     -                 + MW*(eta6*eta7*EL*SBA*
     -                 (-2*Bval(bb00,iintB16(Cha4,Cha5)) + cintB8(Cha5)))
     -                 ))/((-MZ2 + q12)*SW)
               enddo
            enddo
            
            do Neu5 = 1,4
               do Neu4 = 1,4
                  eta1 = ZNeu(Neu5,3)*ZNeuC(Neu4,3) - 
     -                 ZNeu(Neu5,4)*ZNeuC(Neu4,4)
                  eta2 = ZNeu(Neu4,3)*ZNeuC(Neu5,3) - 
     -                 ZNeu(Neu4,4)*ZNeuC(Neu5,4)
                  eta3 = -eta1**2 - eta2**2
                  
                  loopa1 = loopa1 + (-((Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -                 (eta1*eta2*EL*SBA*Bval(bb00,iintB13(Neu4,Neu5)))/
     -                 SW)+(Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(q22*eta1*
     -                 eta2*EL*SBA*Bval(bb1,iintB13(Neu4,Neu5)))/SW +
     -                 (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -                 (eta1*eta2*EL*SBA*cintB7(Neu5))/SW + 
     -                 Alfa/(16.D0*CW2**2*Pi*SW2)*
     -                 (Bval(bb0,iintB13(Neu4,Neu5))*
     -                 (-(MW*(eta3*EL*SBA*MNeu(Neu4)*MNeu(Neu5))) + 
     -                 2*MW*(eta1*eta2*EL*SBA*MNeu2(Neu4))))/SW)/
     -                 (-MZ2 + q22) + (-((Alfa*MW)/(4.D0*CW2**2*Pi*SW2)*
     -                 (eta1*eta2*EL*SBA*
     -                 Bval(bb00,iintB14(Neu4,Neu5)))/SW) + 
     -                 (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*(q12*eta1*eta2*
     -                 EL*SBA*Bval(bb1,iintB14(Neu4,Neu5)))/SW +
     -                 (Alfa*MW)/(8.D0*CW2**2*Pi*SW2)*
     -                 (eta1*eta2*EL*SBA*cintB7(Neu5))/SW + 
     -                 Alfa/(16.D0*CW2**2*Pi*SW2)*
     -                 (Bval(bb0,iintB14(Neu4,Neu5))*
     -                 (-(MW*(eta3*EL*SBA*MNeu(Neu4)*MNeu(Neu5))) + 
     -                 2*MW*(eta1*eta2*EL*SBA*MNeu2(Neu4))))/SW)/
     -                 (-MZ2 + q12)
               enddo
            enddo
         end if
      end if
      

      if (debug .ne. 0) then
         write(*,*)'ZZ loopa1=', loopa1
      end if


** Now we want the form factors themselves         
      a1slfZ = loopa1 + ct 
      
      if (debug .ne. 0) then
         write(*,*)"ZZ slf a1 =", a1slfZ
         write(*,*)"ZZ FERMION contn to slf a1 =", a1selfZfm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the photon-Z self 
*** energy

      subroutine az_selfSUSY(a1slfAZ,q12,q22,q32,q1q2,debug,fm)

      implicit none


** The momenta: q1 and q2 are the momenta of the w bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Dummy variables used when calculating the loop elements
      double complex dup1, dup2
      double complex eta1
      double complex kci1

** The arguments of the integrals (suffix: L -> lepton, B -> boson)
      double complex cint1(3), cint2(3), cint3(3,2), cint4(3,2)
      integer*8 iint1(3), iint2(3), iint3(3,2), iint4(3,2)
      double complex cintL1(3), cintL2(3,2)
      integer*8 iintL1(3), iintL2(3,2)
      double complex cintB1, cintB2, cintB3(2)
      integer*8 iintB1, iintB2, iintB3(2)

** These are the contributions from the counter-term (g-mu-nu), the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1


** This is the formfactor (a1: g-mu-nu) for photon-Z-higgs.  
      double complex a1slfAZ


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for fermions.  'fm' = 1 for AZ and 2 for AZ2
      integer fm
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


** Particle generation do loop variables
      integer Gen4, Sfe4, Sfe5, g4ini
      integer Neu4, Neu5, Cha4, Cha5

** Identity matrix
      double precision IndexDelta(4,4)
      double precision Alfa, Alfa2

      integer debug

** Model parameters
#include "mssm.inc" 
#include "global.inc" 

** This is needed for looptools. 
#include "looptools.h"


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions, 3 or 4 for all SM-type
** particles and 5 or 6 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else 
         g4ini = 1
      end if


** Setting the identity matrix
      do Neu4 = 1, 4
         do Neu5 = 1, 4
            if (Neu4 .eq. Neu5) then
               IndexDelta(Neu4,Neu5) = 1D0
            else
               IndexDelta(Neu4,Neu5) = 0D0
            end if
         end do
      end do

      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Now for the arguments of the integrals.  Note that these are the same
** for both the heavy and light higgs contributions.
** First the quarks and squarks
      do Gen4 = g4ini, 3  
         cint1(Gen4) = A0(Mf2(3,Gen4))
         cint2(Gen4) = A0(Mf2(4,Gen4))
         iint1(Gen4) = Bget(q12,Mf2(3,Gen4),Mf2(3,Gen4))
         iint2(Gen4) = Bget(q12,Mf2(4,Gen4),Mf2(4,Gen4))            
      enddo

      do Sfe4 = 1,2
         do Gen4 = g4ini, 3
            cint3(Gen4,Sfe4) = A0(MSf2(Sfe4,3,Gen4))
            cint4(Gen4,Sfe4) = A0(MSf2(Sfe4,4,Gen4))
            iint3(Gen4,Sfe4) = 
     -           Bget(q12,MSf2(Sfe4,3,Gen4),MSf2(Sfe4,3,Gen4))
            iint4(Gen4,Sfe4) = 
     -           Bget(q12,MSf2(Sfe4,4,Gen4),MSf2(Sfe4,4,Gen4))
         enddo
      enddo

** Now the leptons and sleptons
      if (sector .ne. 1) then
         do Gen4 = 1,3
            cintL1(Gen4) = A0(Mf2(2,Gen4))
            iintL1(Gen4) = Bget(q12,Mf2(2,Gen4),Mf2(2,Gen4))
         enddo

         do Sfe4 = 1,2
            do Gen4 = 1,3
               cintL2(Gen4,Sfe4) = A0(MSf2(Sfe4,2,Gen4))
               iintL2(Gen4,Sfe4) = 
     -              Bget(q12,MSf2(Sfe4,2,Gen4),MSf2(Sfe4,2,Gen4))
            enddo
         enddo
      end if

** Integrals for SM bosons and their superpartners ...
      if (sector .ge. 3) then
         cintB1 = A0(MHp2)
         cintB2 = A0(MW2)
         iintB1 = Bget(q12,MW2,MW2)
         iintB2 = Bget(q12,MHp2,MHp2)
         
** neutralino / chargino integrals
         if (sector .ge. 5) then
            do Cha4 = 1,2
               cintB3(Cha4) = A0(MCha2(Cha4))
               iintB3(Cha4) = Bget(q12,MCha2(Cha4),MCha2(Cha4))
            enddo
         end if
         
      end if


** Clearing the values of each contribution
      loopa1 = (0D0,0D0)  
      ct = (0D0,0D0)   


** Now the insanely long and hideously complicated form factors.  In bits.  

** First the counterterm 
      ct = 1/(2.D0*CW2)*(q12*(-(MW*(dZAZ1*EL*SBA)) -MW*(dZZA1*EL*SBA)) + 
     -     MW*MZ2*(dZZA1*EL*SBA))/((-MZ2 + q12)*SW)

** Fermion contn to the counterterm
      ctFM = 1/(2.D0*CW2)*(q12*(-(MW*(dZAZ1fm*EL*SBA)) -MW*(dZZA1fm*EL*SBA)) + 
     -     MW*MZ2*(dZZA1fm*EL*SBA))/((-MZ2 + q12)*SW)

      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'photon-Z self energy ct=', ct
         write(*,*)'fermion contn to photon-Z self energy ct=', ctFM
      end if
         

** Now the loop: first the quarks and squarks
      do Gen4 = g4ini, 3
         loopa1 = loopa1 + 1/(12.D0*CW*CW2*Pi*SW2)*(Alfa*MW*(q12*EL*SBA*
     -        (6*Bval(bb1,iint1(Gen4)) + 3*Bval(bb1,iint2(Gen4)) + SW2*
     -        (-16*Bval(bb1,iint1(Gen4)) - 4*Bval(bb1,iint2(Gen4))))) + 
     -        Alfa*MW*(EL*SBA*(-12*Bval(bb00,iint1(Gen4)) - 6*
     -        Bval(bb00,iint2(Gen4)) + 6*cint1(Gen4) + 3*cint2(Gen4) + 
     -        SW2*(32*Bval(bb00,iint1(Gen4)) +8*Bval(bb00,iint2(Gen4)) - 
     -        16*cint1(Gen4) - 4*cint2(Gen4)))))/(-MZ2 + q12)
      enddo
      
      do Sfe4 = 1,2
         do Gen4 = g4ini, 3
            dup1 = SW2*(4*USf(Sfe4,1,3,Gen4)*USfC(Sfe4,1,3,Gen4) + 
     -           4*USf(Sfe4,2,3,Gen4)*USfC(Sfe4,2,3,Gen4)) - 
     -           3*USf(Sfe4,1,3,Gen4)*USfC(Sfe4,1,3,Gen4)
            dup2 = SW2*(2*USf(Sfe4,1,4,Gen4)*USfC(Sfe4,1,4,Gen4) + 
     -           2*USf(Sfe4,2,4,Gen4)*USfC(Sfe4,2,4,Gen4)) - 
     -           3*USf(Sfe4,1,4,Gen4)*USfC(Sfe4,1,4,Gen4)

            loopa1 = loopa1 + Alfa/(12.D0*CW*CW2*Pi*SW2)*(MW*(dup1*EL*
     -           SBA*(-4*Bval(bb00,iint3(Gen4,Sfe4)) + 2*
     -           cint3(Gen4,Sfe4))) + MW*(dup2*EL*SBA*(-2*
     -           Bval(bb00,iint4(Gen4,Sfe4)) + cint4(Gen4,Sfe4))))/
     -           (-MZ2 + q12)
         enddo
      enddo
            

** Now the leptons and sleptons
      if (sector .ne. 1) then
            do Gen4 = 1,3
               loopa1 = loopa1 + 1/(4.D0*CW*CW2*Pi*SW2)*(Alfa*MW*(q12*
     -              EL*SBA*(Bval(bb1,iintL1(Gen4)) - 4*SW2*
     -              Bval(bb1,iintL1(Gen4)))) + Alfa*MW*(EL*SBA*(-2*
     -              Bval(bb00,iintL1(Gen4)) + cintL1(Gen4) + SW2*(8*
     -              Bval(bb00,iintL1(Gen4)) - 4*cintL1(Gen4)))))/(-MZ2 + 
     -              q12)
            enddo
            
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  eta1 = SW2*(2*USf(Sfe4,1,2,Gen4)*USfC(Sfe4,1,2,Gen4) + 
     -                 2*USf(Sfe4,2,2,Gen4)*USfC(Sfe4,2,2,Gen4)) - 
     -                 USf(Sfe4,1,2,Gen4)*USfC(Sfe4,1,2,Gen4)

                  loopa1 = loopa1 + (Alfa*MW)/(4.D0*CW*CW2*Pi*SW2)*
     -                 (eta1*EL*SBA*(-2*Bval(bb00,iintL2(Gen4,Sfe4)) + 
     -                 cintL2(Gen4,Sfe4)))/(-MZ2 + q12)
               enddo
            enddo
      end if


** The fermion contn to the form factors 
      if (fm .eq. 1) then
         a1selfAZfm = loopa1 + ctFM
      else if (fm .eq. 2) then
         a1selfAZ2fm = loopa1 + ctFM
      end if


** And the sm bosons and their superpartners ...
      if (sector .ge. 3) then         
         loopa1 = loopa1 + (Bval(bb0,iintB1)*((5*Alfa*CW*MW)/(4.D0*CW2*
     -        Pi*SW2)*(q12*EL*SBA) + (Alfa*MW*MW2)/(2.D0*CW*CW2*Pi*SW2)*
     -        (EL*SBA)) + (Alfa*MW*(-CW2 + SW2))/(4.D0*CW*CW2*Pi*SW2)*
     -        (cintB1*EL*SBA) - (Alfa*MW*(-CW2 + SW2))/(2.D0*CW*CW2*Pi*
     -        SW2)*(EL*SBA*Bval(bb00,iintB2)) + (Alfa*CW*MW)/(2.D0*CW2*
     -        Pi*SW2)*(q12*EL*SBA*Bval(bb1,iintB1)) + Alfa/(2.D0*CW*CW2*
     -        Pi*SW2)*(Bval(bb00,iintB1)*(5*CW2*MW*(EL*SBA) - MW*SW2*
     -        (EL*SBA)))+Alfa/(4.D0*CW*CW2*Pi*SW2)*(cintB2*(-(5*CW2*MW*
     -        (EL*SBA)) + MW*SW2*(EL*SBA))))/(-MZ2 + q12)
         
* charginos and neutralinos
         if (sector .ge. 5) then
            do Cha4 = 1,2
               kci1 = 4*SW2 - 2*UCha(Cha4,1)*UChaC(Cha4,1) - 
     -              UCha(Cha4,2)*UChaC(Cha4,2) - 2*VCha(Cha4,1)*
     -              VChaC(Cha4,1) - VCha(Cha4,2)*VChaC(Cha4,2)

               loopa1 = loopa1 + 1/(4.D0*CW*CW2*Pi*SW2)*(-(Alfa*MW*(q12*
     -              kci1*EL*SBA*Bval(bb1,iintB3(Cha4)))) + Alfa*MW*
     -              (kci1*EL*SBA*(2*Bval(bb00,iintB3(Cha4)) - 
     -              cintB3(Cha4))))/(-MZ2 + q12)
            enddo
         end if
      end if
      

      if (debug .ne. 0) then
         write(*,*)'photon-Z self energy loopa1=', loopa1
      end if


** Now we want the form factors themselves         
      a1slfAZ = loopa1 + ct 
      
      if (debug .ne. 0) then
         write(*,*)"photon-Z slf a1 =", a1slfAZ
         write(*,*)"FERMION contn to photon-Z slf a1 =", a1selfAZfm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the Z -> Higgs
*** mixing in SUSY

      subroutine zh_prpSUSY(a1zh,q12,q22,q32,q1q2,debug,higgno)

      implicit none


** The momenta: q1 and q2 are the momenta of the z bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Integral arguments (suffix L -> 'leptonic', B -> 'bosonic')
      integer*8 iint1(3), iint2(3), iint3(3,2,2), iint4(3,2,2)
      integer*8 iintL1(3), iintL2(3), iintL3(3,2,2)
      integer*8 iintB1, iintB2, iintB3(4,4), iintB4(2,2), iintB5, iintB6

** These are the contributions from the counter-term (g-mu-nu) and the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for WW-higgs.  The first entry in each array is the light higgs, the second 
** entry is the heavy higgs.
      double complex a1zh
      double complex a1zh0, a1zHH, a1zA0


** Particle do loop parameters
      integer Gen4, Sfe4, Sfe5, Cha4, Cha5, Neu4, Neu5, g4ini

** Dummy parameters: dup for (s)quarks, lup for (s)leptons and bup for SM-type
** bosons and their superpartners
      double complex dup1, dup2, dup3, dup4, dup5, dup6, dup7
      double complex dup8, dup9, dup10, dup11, dup12
      double complex lup1, lup2, lup3, lup4, lup5, lup6
      double complex bup1, bup2, bup3, bup4, bup5, bup6, bup7, bup8
      double complex bup9, bup10, bup11, bup12, bup13, bup14, bup15


** This parameter defines which of the neutral higgs bosons is being considered.
** It is set in vbfnlo.dat
      integer higgno

** Parameters: em coupling and identity matrix
      double precision Alfa, Alfa2
      double precision IndexDelta(4,4)

      integer debug

** This is needed for looptools. 
#include "looptools.h"

** Parameter declarations
#include "mssm.inc"
#include "global.inc"


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F using the input EW_APPROX.  It is
** 1 for generation 3 (s)quarks, 2 for all (s)fermions, 3 or 4 for SM-type
** particles and sfermions and 5 or 6 for all particles.
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      do Neu4 = 1, 4
         do Neu5 = 1, 4
            if (Neu4 .eq. Neu5) then
               IndexDelta(Neu4,Neu5) = 1D0
            else
               IndexDelta(Neu4,Neu5) = 0D0
            end if
         end do
      end do


** Defining the integrals: firstly, quark and squark integrals
      do Gen4 = g4ini, 3
         iint1(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
         iint2(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
      enddo

      do Sfe5 = 1,2
         do Sfe4 = 1,2
            do Gen4 = g4ini,3
               iint3(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,3,Gen4),MSf2(Sfe5,3,Gen4))
               iint4(Gen4,Sfe4,Sfe5) = 
     -              Bget(q32,MSf2(Sfe4,4,Gen4),MSf2(Sfe5,4,Gen4))
            enddo
         enddo
      enddo

** Now integrals for the (s)leptons
      if (sector .ne. 1) then
	do Gen4 = 1,3
           iintL1(Gen4) = Bget(q32,MSf2(1,1,Gen4),MSf2(1,1,Gen4))
           iintL2(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
	enddo

	do Sfe5 = 1,2
           do Sfe4 = 1,2
              do Gen4 = 1,3
                 iintL3(Gen4,Sfe4,Sfe5) = 
     -                Bget(q32,MSf2(Sfe4,2,Gen4),MSf2(Sfe5,2,Gen4))
              enddo
           enddo
	enddo
      end if

** And finally the 'bosonic' integrals (i.e. loops involving the
** SM bosons and their superpartners) ..
      if (sector .ge. 3) then
         iintB1 = Bget(q32,MHp2,MHp2)
         iintB2 = Bget(q32,MW2,MW2)

** 'bosonic' integrals used only for production of A0
         if (higgno .eq. 3) then
            iintB1 = Bget(q32,MA02,Mh02)
            iintB2 = Bget(q32,MA02,MHH2)
            iintB5 = Bget(q32,Mh02,MZ2)
            iintB6 = Bget(q32,MHH2,MZ2)
         end if

** chargino and neutralino integrals
         if (sector .ge. 5) then
            do Cha5 = 1,2
               do Cha4 = 1,2
                  iintB4(Cha4,Cha5) = Bget(q32,MCha2(Cha4),MCha2(Cha5))
               enddo
            enddo
            do Neu5 = 1,4
               do Neu4 = 1,4
                  iintB3(Neu4,Neu5) = Bget(q32,MNeu2(Neu4),MNeu2(Neu5))
               enddo
            enddo
         end if
      end if


** Firstly, the lightest higgs production
      if (higgno .eq. 1) then

** Clearing the loop values
         loopa1 = (0D0,0D0)
         ct = (0D0,0D0)

** Firstly, contributions from quarks and squarks
         do Gen4 = g4ini, 3
            loopa1 = loopa1 + Alfa/(8.D0*MW*Pi*SW2)*
     -           ((q12-q22)*(CA*CB*EL*
     -           (SW2*(8*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) + 
     -           16*Bval(bb1,iint1(Gen4))*Mf2(3,Gen4)) - 
     -           3*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) - 
     -           6*Bval(bb1,iint1(Gen4))*Mf2(3,Gen4)) + 
     -           EL*SA*SB*(SW2*
     -           (4*Bval(bb0,iint2(Gen4))*Mf2(4,Gen4) + 
     -           8*Bval(bb1,iint2(Gen4))*Mf2(4,Gen4)) - 
     -           3*Bval(bb0,iint2(Gen4))*Mf2(4,Gen4) - 
     -           6*Bval(bb1,iint2(Gen4))*Mf2(4,Gen4))))/
     -           (CB*(-MZ2 + q32)*SB*SW)
         enddo

         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = g4ini, 3
                  
                  dup1 = -(2*MW*MZ*SW2*(SAB*SB)) + 3*CW*(CA*Mf2(3,Gen4))
                  dup2 = MW*MZ*SW2*(CB*SAB) - 3*CW*(SA*Mf2(4,Gen4))
                  dup3 = MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -                 6*CW*(CA*Mf2(3,Gen4))
                  dup4 = MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -                 6*CW*(SA*Mf2(4,Gen4))
                  dup5 = SW2*(4*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4) + 
     -                 4*USf(Sfe4,2,3,Gen4)*USfC(Sfe5,2,3,Gen4)) - 
     -                 3*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)
                  dup6 = SW2*(2*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4) + 
     -                 2*USf(Sfe4,2,4,Gen4)*USfC(Sfe5,2,4,Gen4)) - 
     -                 3*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)
                  dup7 = 3*CW*((MUEC*SA + CA*Af(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe4,1,3,Gen4)) + 2*dup1*
     -                 USfC(Sfe4,2,3,Gen4)
                  dup8 = 3*CW*((MUE*SA + CA*AfC(3,Gen4,Gen4))*
     -                 Mf(3,Gen4)*USfC(Sfe4,2,3,Gen4)) + dup3*
     -                 USfC(Sfe4,1,3,Gen4)
                  dup9 = 3*CW*((CA*MUEC + SA*Af(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe4,1,4,Gen4)) - 2*dup2*
     -                 USfC(Sfe4,2,4,Gen4)
                  dup10 = 3*CW*((CA*MUE + SA*AfC(4,Gen4,Gen4))*
     -                 Mf(4,Gen4)*USfC(Sfe4,2,4,Gen4)) + dup4*
     -                 USfC(Sfe4,1,4,Gen4)

                  loopa1 = loopa1 + (-(Alfa/(48.D0*CW*MW*Pi*SW2)*((q12-
     -                 q22)*dup5*EL*Bval(bb0,iint3(Gen4,Sfe4,Sfe5))*
     -                 (dup8*USf(Sfe5,1,3,Gen4) + dup7*
     -                 USf(Sfe5,2,3,Gen4)))/(SB*SW)) - Alfa/(24.D0*CW*
     -                 MW*Pi*SW2)*((q12-q22)*dup5*EL*
     -                 Bval(bb1,iint3(Gen4,Sfe4,Sfe5))*(dup8*
     -                 USf(Sfe5,1,3,Gen4) + dup7*USf(Sfe5,2,3,Gen4)))/
     -                 (SB*SW) - Alfa/(48.D0*CW*MW*Pi*SW2)*((q12-q22)*
     -                 dup6*EL*Bval(bb0,iint4(Gen4,Sfe4,Sfe5))*(dup10*
     -                 USf(Sfe5,1,4,Gen4) + dup9*USf(Sfe5,2,4,Gen4)))/
     -                 (CB*SW) - Alfa/(24.D0*CW*MW*Pi*SW2)*((q12-q22)*
     -                 dup6*EL*Bval(bb1,iint4(Gen4,Sfe4,Sfe5))*(dup10*
     -                 USf(Sfe5,1,4,Gen4) + dup9*USf(Sfe5,2,4,Gen4)))/
     -                 (CB*SW))/(-MZ2 + q32)
               enddo
            enddo
         enddo


** Now the contributions from the leptons and sleptons
         if (sector .ne. 1) then
            do Gen4 = 1,3
               loopa1 = loopa1 + Alfa/(16.D0*CW*MW*Pi*SW2)*((q12-q22)*
     -              (MW*MZ*(CB*EL*SAB*(-Bval(bb0,iintL1(Gen4)) - 2*
     -              Bval(bb1,iintL1(Gen4)))) + CW*(EL*SA*(SW2*
     -              (8*Bval(bb0,iintL2(Gen4))*Mf2(2,Gen4) + 
     -              16*Bval(bb1,iintL2(Gen4))*Mf2(2,Gen4)) - 
     -              2*Bval(bb0,iintL2(Gen4))*Mf2(2,Gen4) - 
     -              4*Bval(bb1,iintL2(Gen4))*Mf2(2,Gen4)))))/
     -              (CB*(-MZ2 + q32)*SW)
            enddo

            do Sfe5 = 1,2
               do Sfe4 = 1,2
                  do Gen4 = 1,3
                     lup1 = -(MW*MZ*SW2*(CB*SAB)) + CW*(SA*Mf2(2,Gen4))
                     lup2 = MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 2*CW*
     -                    (SA*Mf2(2,Gen4))
                     lup3 = SW2*(2*USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4) + 2*USf(Sfe4,2,2,Gen4)*
     -                    USfC(Sfe5,2,2,Gen4)) - USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4)
                     lup4 = CW*((CA*MUEC + SA*Af(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 2*lup1*
     -                    USfC(Sfe4,2,2,Gen4)
                     lup5 = CW*((CA*MUE + SA*AfC(2,Gen4,Gen4))*
     -                    Mf(2,Gen4)*USfC(Sfe4,2,2,Gen4)) + lup2*
     -                    USfC(Sfe4,1,2,Gen4)

                     loopa1 = loopa1 + (-(Alfa/(16.D0*CW*MW*Pi*SW2)*
     -                    ((q12-q22)*lup3*EL*
     -                    Bval(bb0,iintL3(Gen4,Sfe4,Sfe5))*(lup5*
     -                    USf(Sfe5,1,2,Gen4) +lup4*USf(Sfe5,2,2,Gen4)))/
     -                    (CB*SW))- Alfa/(8.D0*CW*MW*Pi*SW2)*((q12-q22)*
     -                    lup3*EL*Bval(bb1,iintL3(Gen4,Sfe4,Sfe5))*
     -                    (lup5*USf(Sfe5,1,2,Gen4) + lup4*
     -                    USf(Sfe5,2,2,Gen4)))/(CB*SW))/(-MZ2 + q32)
                  enddo
               enddo
            enddo
         end if

** And finally the SM bosons and their superpartners
         if (sector .ge. 3) then
            loopa1 = loopa1+(-((Alfa*MW*(CW2 - SW2))/(16.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*EL*Bval(bb0,iintB1)*(C2B*SAB + 2*CW2*SBA))/
     -           SW) - (Alfa*MW*(CW2 - SW2))/(8.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*EL*Bval(bb1,iintB1)*(C2B*SAB + 2*CW2*SBA))/
     -           SW+Alfa/(16.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*Bval(bb0,iintB2)*
     -           (C2B*(CW2*MW*(EL*SAB) - MW*SW2*(EL*SAB)) - 
     -           12*CW2**2*MW*(EL*SBA)))/SW + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*Bval(bb1,iintB2)*
     -           (C2B*(CW2*MW*(EL*SAB) - MW*SW2*(EL*SAB)) - 
     -           12*CW2**2*MW*(EL*SBA)))/SW)/(-MZ2 + q32)

** chargino / neutralino contribution
            if (sector .ge. 5) then
               do Cha5 = 1,2
                  do Cha4 = 1,2
                     bup11 = SW2*IndexDelta(Cha4,Cha5) - 
     -                    1/2.D0*(UCha(Cha4,2)*UChaC(Cha5,2)) - 
     -                    UCha(Cha4,1)*UChaC(Cha5,1)
                     bup12 = SA*UCha(Cha5,2)*VCha(Cha4,1) - 
     -                    CA*UCha(Cha5,1)*VCha(Cha4,2)
                     bup13 = SW2*IndexDelta(Cha4,Cha5) - 
     -                    1/2.D0*(VCha(Cha5,2)*VChaC(Cha4,2)) - 
     -                    VCha(Cha5,1)*VChaC(Cha4,1)
                     bup14 = SA*UChaC(Cha4,2)*VChaC(Cha5,1) - 
     -                    CA*UChaC(Cha4,1)*VChaC(Cha5,2)
                     
                     loopa1 = loopa1 + (Alfa/(2.D0*Pi*SW2)*((q12-q22)*
     -                    (bup11*bup12 + bup13*bup14)*EL*
     -                    Bval(bb0,iintB4(Cha4,Cha5))*MCha(Cha4))/(SW*
     -                    Sqrt(2.D0)) + Alfa/(2.D0*Pi*SW2)*((q12-q22)*
     -                    Bval(bb1,iintB4(Cha4,Cha5))*(bup13*bup14*EL*
     -                    MCha(Cha4) + bup12*bup13*EL*MCha(Cha5) + bup11*
     -                    (bup12*EL*MCha(Cha4) + bup14*EL*MCha(Cha5))))/
     -                    (SW*Sqrt(2.D0)))/(-MZ2 + q32)
                  enddo
               enddo
               
               do Neu5 = 1,4
                  do Neu4 = 1,4
                     bup1 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                     bup2 = CA*(CW*ZNeu(Neu5,2) - SW*ZNeu(Neu5,1))
                     bup3 = CW*(SA*ZNeu(Neu5,2)) - SA*SW*ZNeu(Neu5,1)
                     bup4 = SA*ZNeu(Neu5,3) + CA*ZNeu(Neu5,4)
                     bup5 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                     bup6 = ZNeu(Neu5,3)*ZNeuC(Neu4,3) - ZNeu(Neu5,4)*
     -                    ZNeuC(Neu4,4)
                     bup7 = -(CW*ZNeuC(Neu5,2)) + SW*ZNeuC(Neu5,1)
                     bup8 = SA*ZNeuC(Neu5,3) + CA*ZNeuC(Neu5,4)
                     bup9 = -(ZNeu(Neu4,3)*ZNeuC(Neu5,3)) + 
     -                    ZNeu(Neu4,4)*ZNeuC(Neu5,4)
                     bup10 = bup9*(-(bup1*bup4) + bup3*ZNeu(Neu4,3) + 
     -                    bup2*ZNeu(Neu4,4)) - 
     -                    bup6*(bup5*bup8 + bup7*SA*ZNeuC(Neu4,3) + 
     -                    CA*bup7*ZNeuC(Neu4,4))
                     
                     loopa1 = loopa1 + (-(Alfa/(16.D0*CW*Pi*SW2)*
     -                    ((q12-q22)*bup10*EL*Bval(bb0,iintB3(Neu4,Neu5))*
     -                    MNeu(Neu4))/SW) + Bval(bb1,iintB3(Neu4,Neu5))*
     -                    (-(Alfa/(16.D0*CW*Pi*SW2)*
     -                    ((q12-q22)*bup10*EL*MNeu(Neu4))/SW) + 
     -                    Alfa/(16.D0*CW*Pi*SW2)*
     -                    ((q12-q22)*EL*MNeu(Neu5)*(-(bup6*
     -                    (-(bup1*bup4) + bup3*ZNeu(Neu4,3) + 
     -                    bup2*ZNeu(Neu4,4))) + 
     -                    bup9*(bup5*bup8 + bup7*SA*ZNeuC(Neu4,3) + 
     -                    CA*bup7*ZNeuC(Neu4,4))))/SW))/
     -                    (-MZ2 + q32)
                  enddo
               enddo
            end if
         end if

** End of lightest higgs production
      end if


** And now for the heavier higgs production
      if (higgno .eq. 2) then
         
** Clearing loopa1
         loopa1 = (0D0,0D0)
         ct = (0D0,0D0)

** First the quarks and squarks
         do Gen4 = g4ini, 3
            loopa1 = loopa1 + Alfa/(8.D0*MW*Pi*SW2)*
     -           ((q12-q22)*(CB*EL*SA*
     -           (SW2*(8*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) + 
     -           16*Bval(bb1,iint1(Gen4))*Mf2(3,Gen4)) - 
     -           3*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) - 
     -           6*Bval(bb1,iint1(Gen4))*Mf2(3,Gen4)) + 
     -           CA*EL*SB*(SW2*
     -           (-4*Bval(bb0,iint2(Gen4))*Mf2(4,Gen4) - 
     -           8*Bval(bb1,iint2(Gen4))*Mf2(4,Gen4)) + 
     -           3*Bval(bb0,iint2(Gen4))*Mf2(4,Gen4) + 
     -           6*Bval(bb1,iint2(Gen4))*Mf2(4,Gen4))))/
     -           (CB*(-MZ2 + q32)*SB*SW)
         enddo

         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = g4ini,3
                  dup1 = -(CA*MUEC) + SA*Af(3,Gen4,Gen4)
                  dup2 = -(MUEC*SA) + CA*Af(4,Gen4,Gen4)
                  dup3 = -(CA*MUE) + SA*AfC(3,Gen4,Gen4)
                  dup4 = -(MUE*SA) + CA*AfC(4,Gen4,Gen4)
                  dup5 = 2*MW*MZ*SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen4))
                  dup6 = MW*MZ*SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen4))
                  dup7 = MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -                 6*CW*(SA*Mf2(3,Gen4))
                  dup8 = MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -                 6*CW*(CA*Mf2(4,Gen4))
                  dup9 = SW2*(4*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4) + 
     -                 4*USf(Sfe4,2,3,Gen4)*USfC(Sfe5,2,3,Gen4)) - 
     -                 3*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)
                  dup10 = SW2*(2*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)+ 
     -                 2*USf(Sfe4,2,4,Gen4)*USfC(Sfe5,2,4,Gen4)) - 
     -                 3*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)
                  dup11 = USf(Sfe5,1,3,Gen4)*
     -                 (3*CW*(dup3*Mf(3,Gen4)*USfC(Sfe4,2,3,Gen4)) + 
     -                 dup7*USfC(Sfe4,1,3,Gen4)) + 
     -                 USf(Sfe5,2,3,Gen4)*
     -                 (3*CW*(dup1*Mf(3,Gen4)*USfC(Sfe4,1,3,Gen4)) + 
     -                 2*dup5*USfC(Sfe4,2,3,Gen4))
                  dup12 = USf(Sfe5,1,4,Gen4)*
     -                 (3*CW*(dup4*Mf(4,Gen4)*USfC(Sfe4,2,4,Gen4)) + 
     -                 dup8*USfC(Sfe4,1,4,Gen4)) + 
     -                 USf(Sfe5,2,4,Gen4)*
     -                 (3*CW*(dup2*Mf(4,Gen4)*USfC(Sfe4,1,4,Gen4)) - 
     -                 2*dup6*USfC(Sfe4,2,4,Gen4))

                  loopa1 = loopa1 + Alfa/(48.D0*CW*MW*Pi*SW2)*
     -                 ((q12-q22)*(CB*dup11*dup9*EL*
     -                 (-Bval(bb0,iint3(Gen4,Sfe4,Sfe5)) - 
     -                 2*Bval(bb1,iint3(Gen4,Sfe4,Sfe5))) + 
     -                 dup10*dup12*EL*SB*
     -                 (Bval(bb0,iint4(Gen4,Sfe4,Sfe5)) + 
     -                 2*Bval(bb1,iint4(Gen4,Sfe4,Sfe5)))))/
     -                 (CB*(-MZ2 + q32)*SB*SW)
               enddo
            enddo
         enddo

** Secondly the leptons and sleptons
         if (sector .ne. 1) then
            do Gen4 = 1,3
               loopa1 = loopa1 + Alfa/(16.D0*CW*MW*Pi*SW2)*
     -              ((q12-q22)*(MW*MZ*(CAB*CB*EL*
     -              (Bval(bb0,iintL1(Gen4)) + 2*Bval(bb1,iintL1(Gen4)))
     -              ) + CW*(CA*EL*
     -              (SW2*(-8*Bval(bb0,iintL2(Gen4))*Mf2(2,Gen4) - 
     -              16*Bval(bb1,iintL2(Gen4))*Mf2(2,Gen4)) + 
     -              2*Bval(bb0,iintL2(Gen4))*Mf2(2,Gen4) + 
     -              4*Bval(bb1,iintL2(Gen4))*Mf2(2,Gen4)))))/
     -              (CB*(-MZ2 + q32)*SW)
            enddo

            do Sfe5 = 1,2
               do Sfe4 = 1,2
                  do Gen4 = 1,3
                     lup1 = MUEC*SA - CA*Af(2,Gen4,Gen4)
                     lup2 = MUE*SA - CA*AfC(2,Gen4,Gen4)
                     lup3 = MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen4))
                     lup4 = MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                    2*CW*(CA*Mf2(2,Gen4))
                     lup5 = SW2*(2*USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4) + 
     -                    2*USf(Sfe4,2,2,Gen4)*USfC(Sfe5,2,2,Gen4)) - 
     -                    USf(Sfe4,1,2,Gen4)*USfC(Sfe5,1,2,Gen4)
                     lup6 = USf(Sfe5,1,2,Gen4)*
     -                    (CW*(lup2*Mf(2,Gen4)*USfC(Sfe4,2,2,Gen4)) + 
     -                    lup4*USfC(Sfe4,1,2,Gen4)) + 
     -                    USf(Sfe5,2,2,Gen4)*
     -                    (CW*(lup1*Mf(2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 
     -                    2*lup3*USfC(Sfe4,2,2,Gen4))
                     
                     loopa1 = loopa1 + Alfa/(16.D0*CW*MW*Pi*SW2)*
     -                    ((q12-q22)*lup5*lup6*EL*
     -                    (-Bval(bb0,iintL3(Gen4,Sfe4,Sfe5)) - 
     -                    2*Bval(bb1,iintL3(Gen4,Sfe4,Sfe5))))/
     -                    (CB*(-MZ2 + q32)*SW)
                  enddo
               enddo
            enddo
         end if

** Finally the SM bosons and their superpartners
         if (sector .ge. 3) then
            loopa1 = loopa1 + ((Alfa*MW*(CW2 - SW2))/(16.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*EL*Bval(bb0,iintB1)*(C2B*CAB - 2*CW2*CBA))/
     -           SW+(Alfa*MW*(CW2 - SW2))/(8.D0*CW2*Pi*SW2)*((q12-q22)*
     -           EL*Bval(bb1,iintB1)*(C2B*CAB - 2*CW2*CBA))/SW+
     -           Alfa/(16.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*Bval(bb0,iintB2)*
     -           (C2B*CAB*(-(CW2*MW*EL) + MW*SW2*EL) - 
     -           12*CW2**2*MW*(CBA*EL)))/SW + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*Bval(bb1,iintB2)*
     -           (C2B*CAB*(-(CW2*MW*EL) + MW*SW2*EL) - 
     -           12*CW2**2*MW*(CBA*EL)))/SW)/(-MZ2 + q32)

** chargino / neutralino contribution
            if (sector .ge. 5) then
               do Cha5 = 1,2
                  do Cha4 = 1,2
                     bup12 = SW2*IndexDelta(Cha4,Cha5) - 
     -                    1/2.D0*(UCha(Cha4,2)*UChaC(Cha5,2)) - 
     -                    UCha(Cha4,1)*UChaC(Cha5,1)
                     bup13 = CA*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    SA*UCha(Cha5,1)*VCha(Cha4,2)
                     bup14 = SW2*IndexDelta(Cha4,Cha5) - 
     -                    1/2.D0*(VCha(Cha5,2)*VChaC(Cha4,2)) - 
     -                    VCha(Cha5,1)*VChaC(Cha4,1)
                     bup15 = CA*UChaC(Cha4,2)*VChaC(Cha5,1) + 
     -                    SA*UChaC(Cha4,1)*VChaC(Cha5,2)
                     
                     loopa1 = loopa1 + (-(Alfa/(2.D0*Pi*SW2)*
     -                    ((q12-q22)*(bup12*bup13 + bup14*bup15)*EL*
     -                    Bval(bb0,iintB4(Cha4,Cha5))*MCha(Cha4))/
     -                    (SW*Sqrt(2.D0))) + Alfa/(2.D0*Pi*SW2)*
     -                    ((q12-q22)*Bval(bb1,iintB4(Cha4,Cha5))*
     -                    (-(bup14*bup15*EL*MCha(Cha4)) - 
     -                    bup13*bup14*EL*MCha(Cha5) + 
     -                    bup12*(-(bup13*EL*MCha(Cha4)) - 
     -                    bup15*EL*MCha(Cha5))))/(SW*Sqrt(2.D0)))/
     -                    (-MZ2 + q32)
                  enddo
               enddo
               
               do Neu5 = 1,4
                  do Neu4 = 1,4
                     bup1 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                     bup2 = -(CW*ZNeu(Neu5,2)) + SW*ZNeu(Neu5,1)
                     bup3 = CW*(SA*ZNeu(Neu5,2)) - SA*SW*ZNeu(Neu5,1)
                     bup4 = CA*ZNeu(Neu5,3) - SA*ZNeu(Neu5,4)
                     bup5 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                     bup6 = ZNeu(Neu5,3)*ZNeuC(Neu4,3) - ZNeu(Neu5,4)*
     -                    ZNeuC(Neu4,4)
                     bup7 = -(CW*ZNeuC(Neu5,2)) + SW*ZNeuC(Neu5,1)
                     bup8 = CW*(SA*ZNeuC(Neu5,2)) - SA*SW*ZNeuC(Neu5,1)
                     bup9 = CA*ZNeuC(Neu5,3) - SA*ZNeuC(Neu5,4)
                     bup10 = -(ZNeu(Neu4,3)*ZNeuC(Neu5,3)) + 
     -                    ZNeu(Neu4,4)*ZNeuC(Neu5,4)
                     bup11 = bup10*(bup1*bup4 + CA*bup2*ZNeu(Neu4,3) + 
     -                    bup3*ZNeu(Neu4,4)) + 
     -                    bup6*(bup5*bup9 + CA*bup7*ZNeuC(Neu4,3) + 
     -                    bup8*ZNeuC(Neu4,4))
                     
                     loopa1 = loopa1 + (-(Alfa/(16.D0*CW*Pi*SW2)*
     -                    ((q12-q22)*bup11*EL*Bval(bb0,iintB3(Neu4,Neu5))*
     -                    MNeu(Neu4))/SW) + 
     -                    Bval(bb1,iintB3(Neu4,Neu5))*
     -                    (-(Alfa/(16.D0*CW*Pi*SW2)*
     -                    ((q12-q22)*bup11*EL*MNeu(Neu4))/SW) - 
     -                    Alfa/(16.D0*CW*Pi*SW2)*
     -                    ((q12-q22)*EL*MNeu(Neu5)*
     -                    (bup6*(bup1*bup4 + CA*bup2*ZNeu(Neu4,3) + 
     -                    bup3*ZNeu(Neu4,4)) + 
     -                    bup10*
     -                    (bup5*bup9 + CA*bup7*ZNeuC(Neu4,3) + 
     -                    bup8*ZNeuC(Neu4,4))))/SW))/(-MZ2 + q32)
                  enddo
               enddo
            end if
         end if

** End of second heaviest higgs production
      end if


** Lastly, calculation of A0 formfactors
      if (higgno .eq. 3) then

         ct = (0D0,0D0)
         loopa1 = (0D0,0D0)

** Counterterm
         ct =  -CW*((q12-q22)*dMAZsq1*EL)/((-MZ2 + q32)*SW)

** (S)quark contributions
         do Gen4 = g4ini, 3
            loopa1 = loopa1 + (Alfa*cI)/(8.D0*MW*Pi*SW2)*
     -           ((q12-q22)*EL*(3*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) - 
     -           3*TB2*Bval(bb0,iint2(Gen4))*Mf2(4,Gen4)))/
     -           ((-MZ2 + q32)*SW*TB)
         enddo

         do Sfe5 = 1,2
            do Sfe4 = 1,2
               do Gen4 = g4ini,3
                  dup1 = -((MUEC*TB + Af(3,Gen4,Gen4))*
     -                 USf(Sfe5,2,3,Gen4)*USfC(Sfe4,1,3,Gen4)) + 
     -                 (MUE*TB + AfC(3,Gen4,Gen4))*USf(Sfe5,1,3,Gen4)*
     -                 USfC(Sfe4,2,3,Gen4)
                  dup2 = 4*SW2*(USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4) + 
     -                 USf(Sfe4,2,3,Gen4)*USfC(Sfe5,2,3,Gen4)) - 
     -                 3*USf(Sfe4,1,3,Gen4)*USfC(Sfe5,1,3,Gen4)
                  dup3 = 2*SW2*(USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4) + 
     -                 USf(Sfe4,2,4,Gen4)*USfC(Sfe5,2,4,Gen4)) - 
     -                 3*USf(Sfe4,1,4,Gen4)*USfC(Sfe5,1,4,Gen4)
                  dup4 = -((MUEC + TB*Af(4,Gen4,Gen4))*
     -                 USf(Sfe5,2,4,Gen4)*USfC(Sfe4,1,4,Gen4)) + 
     -                 (MUE + TB*AfC(4,Gen4,Gen4))*USf(Sfe5,1,4,Gen4)*
     -                 USfC(Sfe4,2,4,Gen4)
                  
                  loopa1 = loopa1 + (Alfa*cI)/(16.D0*MW*Pi*SW2)*
     -                 ((q12-q22)*EL*(dup1*dup2*
     -                 (Bval(bb0,iint3(Gen4,Sfe4,Sfe5)) + 
     -                 2*Bval(bb1,iint3(Gen4,Sfe4,Sfe5)))*Mf(3,Gen4) + 
     -                 dup3*dup4*TB*
     -                 (-Bval(bb0,iint4(Gen4,Sfe4,Sfe5)) - 
     -                 2*Bval(bb1,iint4(Gen4,Sfe4,Sfe5)))*Mf(4,Gen4)))/
     -                 ((-MZ2 + q32)*SW*TB)
               enddo
            enddo
         enddo

** (S)lepton contributions
         if (sector .ne. 1) then
            do Gen4 = 1,3
               loopa1 = loopa1 - (Alfa*cI)/(8.D0*MW*Pi*SW2)*((q12-q22)*
     -              EL*TB*Bval(bb0,iintL2(Gen4))*Mf2(2,Gen4))/
     -              ((-MZ2 + q32)*SW)
            enddo
            
            do Sfe5 = 1,2
               do Sfe4 = 1,2
                  do Gen4 = 1,3
                     lup1 = 2*SW2*(USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4) + USf(Sfe4,2,2,Gen4)*
     -                    USfC(Sfe5,2,2,Gen4)) - USf(Sfe4,1,2,Gen4)*
     -                    USfC(Sfe5,1,2,Gen4)
                     
                     lup2 = -((MUEC + TB*Af(2,Gen4,Gen4))*
     -                    USf(Sfe5,2,2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 
     -                    (MUE + TB*AfC(2,Gen4,Gen4))*
     -                    USf(Sfe5,1,2,Gen4)*USfC(Sfe4,2,2,Gen4)
                     
                     loopa1 = loopa1 + (Alfa*cI)/(16.D0*MW*Pi*SW2)*
     -                    ((q12-q22)*lup1*lup2*EL*
     -                    (-Bval(bb0,iintL3(Gen4,Sfe4,Sfe5)) - 
     -                    2*Bval(bb1,iintL3(Gen4,Sfe4,Sfe5)))*
     -                    Mf(2,Gen4))/((-MZ2 + q32)*SW)
                  enddo
               enddo
            enddo
         end if

** Contribution from SM-type bosons and their superpartners
         if (sector .ge. 3) then
            loopa1 = loopa1 + (Alfa*cI*MW)/(16.D0*CW2*Pi*SW2)*
     -           ((q12-q22)*EL*(C2B*
     -           (CBA*SAB*(Bval(bb0,iintB1) + 2*Bval(bb1,iintB1)) + 
     -           CAB*SBA*(Bval(bb0,iintB2) + 2*Bval(bb1,iintB2))) + 
     -           CAB*CBA*S2B*(Bval(bb0,iintB6) + 2*Bval(bb1,iintB6)) + 
     -           SBA*(S2B*SAB*
     -           (-Bval(bb0,iintB5) - 2*Bval(bb1,iintB5)) + 
     -           CBA*(2*(Bval(bb0,iintB6) + Bval(bb1,iintB5)) - 
     -           2*(Bval(bb0,iintB5) + Bval(bb1,iintB6))))))/
     -           ((-MZ2 + q32)*SW)

** chargino / neutralino contribution
            if (sector .ge. 5) then
               do Cha5 = 1,2
                  do Cha4 = 1,2
                     bup12 = SW2*IndexDelta(Cha4,Cha5) - 
     -                    1/2.D0*(UCha(Cha4,2)*UChaC(Cha5,2)) - 
     -                    UCha(Cha4,1)*UChaC(Cha5,1)
                     bup13 = SB*UCha(Cha5,2)*VCha(Cha4,1) + 
     -                    CB*UCha(Cha5,1)*VCha(Cha4,2)
                     bup14 = SW2*IndexDelta(Cha4,Cha5) - 
     -                    1/2.D0*(VCha(Cha5,2)*VChaC(Cha4,2)) - 
     -                    VCha(Cha5,1)*VChaC(Cha4,1)
                     bup15 = SB*UChaC(Cha4,2)*VChaC(Cha5,1) + 
     -                    CB*UChaC(Cha4,1)*VChaC(Cha5,2)
                     
                     loopa1 = loopa1 + (Alfa*cI)/(Pi*SW2)*((q12-q22)*EL*
     -                    (-(1/2.D0*((bup12*bup13 - bup14*bup15)*
     -                    Bval(bb0,iintB4(Cha4,Cha5))*MCha(Cha4))) + 
     -                    1/2.D0*(Bval(bb1,iintB4(Cha4,Cha5))*
     -                    (bup14*(bup15*MCha(Cha4) - bup13*MCha(Cha5)) + 
     -                    bup12*(-(bup13*MCha(Cha4)) + bup15*MCha(Cha5)))
     -                    )))/((-MZ2 + q32)*SW*Sqrt(2.D0))
                  enddo
               enddo
               
               do Neu5 = 1,4
                  do Neu4 = 1,4
                     bup1 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                     bup2 = -(CW*ZNeu(Neu5,2)) + SW*ZNeu(Neu5,1)
                     bup3 = CB*(CW*ZNeu(Neu5,2) - SW*ZNeu(Neu5,1))
                     bup4 = SB*ZNeu(Neu5,3) - CB*ZNeu(Neu5,4)
                     bup5 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                     bup6 = ZNeu(Neu5,3)*ZNeuC(Neu4,3) - ZNeu(Neu5,4)*
     -                    ZNeuC(Neu4,4)
                     bup7 = -(CW*ZNeuC(Neu5,2)) + SW*ZNeuC(Neu5,1)
                     bup8 = CB*(CW*ZNeuC(Neu5,2) - SW*ZNeuC(Neu5,1))
                     bup9 = SB*ZNeuC(Neu5,3) - CB*ZNeuC(Neu5,4)
                     bup10 = -(ZNeu(Neu4,3)*ZNeuC(Neu5,3)) + 
     -                    ZNeu(Neu4,4)*ZNeuC(Neu5,4)
                     bup11 = bup10*(bup1*bup4 + bup2*SB*ZNeu(Neu4,3) + 
     -                    bup3*ZNeu(Neu4,4)) - 
     -                    bup6*(bup5*bup9 + bup7*SB*ZNeuC(Neu4,3) + 
     -                    bup8*ZNeuC(Neu4,4))
                     
                     loopa1 = loopa1 + (Alfa*cI)/(CW*Pi*SW2)*
     -                    ((q12-q22)*EL*(-(1/16.D0*
     -                    (bup11*Bval(bb0,iintB3(Neu4,Neu5))*MNeu(Neu4))) + 
     -                    Bval(bb1,iintB3(Neu4,Neu5))*
     -                    (-(1/16.D0*(bup11*MNeu(Neu4))) + 
     -                    1/16.D0*(MNeu(Neu5)*
     -                    (-(bup6*
     -                    (bup1*bup4 + bup2*SB*ZNeu(Neu4,3) + 
     -                    bup3*ZNeu(Neu4,4))) + 
     -                    bup10*
     -                    (bup5*bup9 + bup7*SB*ZNeuC(Neu4,3) + 
     -                    bup8*ZNeuC(Neu4,4)))))))/
     -                    ((-MZ2 + q32)*SW)
                  enddo
               enddo
            end if
         end if         

** End of A0 formfactor calculation
      end if


      a1zh = loopa1 + ct


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the photon -> Higgs
*** mixing in SUSY

      subroutine gamh_prpSUSY(a1gamh,q12,q22,q32,q1q2,debug,higgno)

      implicit none


** The momenta: q1 and q2 are the momenta of the z bosons, q3 the momentum of
** the higgs boson and qx2 is the square of qx, and q1q2 is the dot product
** q1.q2. 
      double precision q12, q22, q32, q1q2


** Integral arguments (suffix L -> 'leptonic', B -> 'bosonic')
      integer*8 iint1(3), iint2(3), iint3(3,2), iint4(3,2)
      integer*8 iintL1(3), iintL2(3,2)
      integer*8 iintB1, iintB2, iintB3(2)

** These are the contributions from the counter-term (g-mu-nu) and the loop
** (g-mu-nu, q2-mu q1-nu and levi-civita respectively) 
      double complex ct, loopa1


** These are the formfactors (a1: g-mu-nu, a2: q2-mu q1-nu, a3: levi-civita)
** for WW-higgs.  
      double complex a1gamh
      double complex a1zgam0, a1gamHH, a1gamA0


** Particle do loop parameters
      integer Gen4, Sfe4, Sfe5, Cha4, Cha5, Neu4, Neu5, g4ini

** Dummy parameters
      double complex dup1, dup2, dup3, dup4, dup5, dup6, dup7
      double complex dup8, dup9, dup10
      double complex lup1, lup2, lup3, lup4, lup5
      double complex bup1


** This parameter defines which of the neutral higgs bosons is being considered.
** It is set in vbfnlo.dat
      integer higgno


** Necessary parameters: em coupling, complex i and identity matrix
      double precision Alfa, Alfa2
      double precision IndexDelta(4,4)

      integer debug

** This is needed for looptools. 
#include "looptools.h"

** Parameter declarations
#include "mssm.inc"
#include "global.inc"

** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set on line 80 of koppln.F and is 1 for generation
** 3 quarks, 2 for all fermions, 3 or 4 for SM-type particles plus sfermions and
** 5 or 6 for all particles
      if (sector .eq. 1) then
         g4ini = 3
      else
         g4ini = 1
      end if


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      do Neu4 = 1, 4
         do Neu5 = 1, 4
            if (Neu4 .eq. Neu5) then
               IndexDelta(Neu4,Neu5) = 1D0
            else
               IndexDelta(Neu4,Neu5) = 0D0
            end if
         end do
      end do


** Defining the integrals: firstly, quark and squark integrals
      do Gen4 = 1,3
         iint1(Gen4) = Bget(q32,Mf2(3,Gen4),Mf2(3,Gen4))
         iint2(Gen4) = Bget(q32,Mf2(4,Gen4),Mf2(4,Gen4))
      enddo

      do Sfe4 = 1,2
         do Gen4 = 1,3
            iint3(Gen4,Sfe4) = 
     -           Bget(q32,MSf2(Sfe4,3,Gen4),MSf2(Sfe4,3,Gen4))
            iint4(Gen4,Sfe4) = 
     -           Bget(q32,MSf2(Sfe4,4,Gen4),MSf2(Sfe4,4,Gen4))
         enddo
      enddo

** Now the (s)leptonic integrals
      if (sector .ne. 1) then
         do Gen4 = 1,3
            iintL1(Gen4) = Bget(q32,Mf2(2,Gen4),Mf2(2,Gen4))
         enddo
         do Sfe4 = 1,2
            do Gen4 = 1,3
               iintL2(Gen4,Sfe4) = 
     -              Bget(q32,MSf2(Sfe4,2,Gen4),MSf2(Sfe4,2,Gen4))
            enddo
         enddo
      end if

** And finally the 'bosonic' integrals: i.e. SM-type bosons and their 
** superpartners
      if (sector .ge. 3) then
         iintB1 = Bget(q32,MHp2,MHp2)
         iintB2 = Bget(q32,MW2,MW2)
** neutralino / chargino integrals
         if (sector .ge. 5) then
            do Cha4 = 1,2
               iintB3(Cha4) = Bget(q32,MCha2(Cha4),MCha2(Cha4))
            enddo
         end if
      end if


** Firstly, the formfactor for the light Higgs
      if (higgno .eq. 1) then

** Firstly, contributions from quarks and squarks
         do Gen4 = g4ini,3
            loopa1 = loopa1 + Alfa/(2.D0*MW*Pi)*((q12-q22)*
     -           (CA*CB*EL*(-2*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) - 
     -           4*Bval(bb1,iint1(Gen4))*Mf2(3,Gen4)) + 
     -           EL*SA*SB*(-(Bval(bb0,iint2(Gen4))*Mf2(4,Gen4)) - 
     -           2*Bval(bb1,iint2(Gen4))*Mf2(4,Gen4))))/
     -           (CB*q32*SB*SW)
         enddo

         do Sfe4 = 1,2
            do Gen4 = g4ini, 3
               dup1 = -(2*MW*MZ*SW2*(SAB*SB)) + 3*CW*(CA*Mf2(3,Gen4))
               dup2 = MW*MZ*SW2*(CB*SAB) - 3*CW*(SA*Mf2(4,Gen4))
               dup3 = MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -              6*CW*(CA*Mf2(3,Gen4))
               dup4 = MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -              6*CW*(SA*Mf2(4,Gen4))
               dup5 = 3*CW*((MUEC*SA + CA*Af(3,Gen4,Gen4))*Mf(3,Gen4)*
     -              USfC(Sfe4,1,3,Gen4)) + 2*dup1*USfC(Sfe4,2,3,Gen4)
               dup6 = 3*CW*((MUE*SA + CA*AfC(3,Gen4,Gen4))*Mf(3,Gen4)*
     -              USfC(Sfe4,2,3,Gen4)) + dup3*USfC(Sfe4,1,3,Gen4)
               dup7 = 3*CW*((CA*MUEC + SA*Af(4,Gen4,Gen4))*Mf(4,Gen4)*
     -              USfC(Sfe4,1,4,Gen4)) - 2*dup2*USfC(Sfe4,2,4,Gen4)
               dup8 = 3*CW*((CA*MUE + SA*AfC(4,Gen4,Gen4))*Mf(4,Gen4)*
     -              USfC(Sfe4,2,4,Gen4)) + dup4*USfC(Sfe4,1,4,Gen4)

               loopa1 = loopa1 + 
     -              (Alfa/(12.D0*CW*MW*Pi)*
     -              ((q12-q22)*EL*Bval(bb0,iint3(Gen4,Sfe4))*
     -              (dup6*USf(Sfe4,1,3,Gen4) + dup5*USf(Sfe4,2,3,Gen4))
     -              )/(SB*SW) + 
     -              Alfa/(6.D0*CW*MW*Pi)*
     -              ((q12-q22)*EL*Bval(bb1,iint3(Gen4,Sfe4))*
     -              (dup6*USf(Sfe4,1,3,Gen4) + dup5*USf(Sfe4,2,3,Gen4))
     -              )/(SB*SW) + 
     -              Alfa/(24.D0*CW*MW*Pi)*
     -              ((q12-q22)*EL*Bval(bb0,iint4(Gen4,Sfe4))*
     -              (dup8*USf(Sfe4,1,4,Gen4) + dup7*USf(Sfe4,2,4,Gen4))
     -              )/(CB*SW) + 
     -              Alfa/(12.D0*CW*MW*Pi)*
     -              ((q12-q22)*EL*Bval(bb1,iint4(Gen4,Sfe4))*
     -              (dup8*USf(Sfe4,1,4,Gen4) + dup7*USf(Sfe4,2,4,Gen4))
     -              )/(CB*SW))/q32
            enddo
         enddo


** Now the leptonic contributions         
         if (sector .ne. 1) then
            do Gen4 = 1,3
               loopa1 = loopa1 + 
     -              Alfa/(2.D0*MW*Pi)*((q12-q22)*EL*SA*
     -              (-(Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4)) - 
     -              2*Bval(bb1,iintL1(Gen4))*Mf2(2,Gen4)))/(CB*q32*SW)
            enddo

            do Sfe4 = 1,2
               do Gen4 = 1,3
                  lup1 = -(MW*MZ*SW2*(CB*SAB)) + CW*(SA*Mf2(2,Gen4))
                  lup2 = MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 
     -                 2*CW*(SA*Mf2(2,Gen4))
                  lup3 = CW*((CA*MUEC + SA*Af(2,Gen4,Gen4))*Mf(2,Gen4)*
     -                 USfC(Sfe4,1,2,Gen4)) + 2*lup1*USfC(Sfe4,2,2,Gen4)
                  lup4 = CW*((CA*MUE + SA*AfC(2,Gen4,Gen4))*Mf(2,Gen4)*
     -                 USfC(Sfe4,2,2,Gen4)) + lup2*USfC(Sfe4,1,2,Gen4)

                  loopa1 = loopa1 + 
     -                 (Alfa/(8.D0*CW*MW*Pi)*
     -                 ((q12-q22)*EL*Bval(bb0,iintL2(Gen4,Sfe4))*
     -                 (lup4*USf(Sfe4,1,2,Gen4) + lup3*USf(Sfe4,2,2,Gen4))
     -                 )/(CB*SW) + 
     -                 Alfa/(4.D0*CW*MW*Pi)*
     -                 ((q12-q22)*EL*Bval(bb1,iintL2(Gen4,Sfe4))*
     -                 (lup4*USf(Sfe4,1,2,Gen4) + lup3*USf(Sfe4,2,2,Gen4))
     -                 )/(CB*SW))/q32
               enddo
            enddo

         end if


** Finally, the contributions from the SM-type bosons and their superpartners
         if (sector .ge. 3) then
            loopa1 = loopa1 + (-((Alfa*MW)/(8.D0*CW2*Pi)*((q12-q22)*
     -           EL*Bval(bb0,iintB1)*(C2B*SAB + 2*CW2*SBA))/SW) - 
     -           (Alfa*MW)/(4.D0*CW2*Pi)*((q12-q22)*EL*
     -           Bval(bb1,iintB1)*(C2B*SAB + 2*CW2*SBA))/SW-(Alfa*
     -           MW)/(8.D0*CW2*Pi)*((q12-q22)*EL*Bval(bb0,iintB2)*
     -           (-(C2B*SAB) + 6*CW2*SBA))/SW - (Alfa*MW)/(4.D0*CW2*
     -           Pi)*((q12-q22)*EL*Bval(bb1,iintB2)*(-(C2B*SAB) + 6*
     -           CW2*SBA))/SW)/q32

** neutralino / chargino contribution
            if (sector .ge. 5) then
               do Cha4 = 1,2
                  bup1 = SA*(UCha(Cha4,2)*VCha(Cha4,1) + 
     -                 UChaC(Cha4,2)*VChaC(Cha4,1)) + 
     -                 CA*(-(UCha(Cha4,1)*VCha(Cha4,2)) - 
     -                 UChaC(Cha4,1)*VChaC(Cha4,2))
                  
                  loopa1 = loopa1 + Alfa/(2.D0*Pi)*((q12-q22)*bup1*EL*
     -                 (-(Bval(bb0,iintB3(Cha4))*MCha(Cha4)) - 
     -                 2*Bval(bb1,iintB3(Cha4))*MCha(Cha4)))/
     -                 (q32*SW*Sqrt(2.D0))
               enddo
            end if
         end if

** End of light higgs formfactor calculation
      end if

      
** Calculation of the second heaviest Higgs calculation h2
      if (higgno .eq. 2) then

** Starting with contributions for (s)quarks
         do Gen4 = g4ini, 3
            loopa1 = loopa1 + Alfa/(2.D0*MW*Pi)*((q12-q22)*
     -           (CB*EL*SA*(-2*Bval(bb0,iint1(Gen4))*Mf2(3,Gen4) - 
     -           4*Bval(bb1,iint1(Gen4))*Mf2(3,Gen4)) + 
     -           CA*EL*SB*(Bval(bb0,iint2(Gen4))*Mf2(4,Gen4) + 
     -           2*Bval(bb1,iint2(Gen4))*Mf2(4,Gen4))))/
     -           (CB*q32*SB*SW)
         enddo

         do Sfe4 = 1,2
            do Gen4 = g4ini, 3

               dup1 = -(CA*MUEC) + SA*Af(3,Gen4,Gen4)
               dup2 = -(MUEC*SA) + CA*Af(4,Gen4,Gen4)
               dup3 = -(CA*MUE) + SA*AfC(3,Gen4,Gen4)
               dup4 = -(MUE*SA) + CA*AfC(4,Gen4,Gen4)
               dup5 = 2*MW*MZ*SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen4))
               dup6 = MW*MZ*SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen4))
               dup7 = MW*MZ*(3 - 4*SW2)*(CAB*SB) + 6*CW*(SA*Mf2(3,Gen4))
               dup8 = MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -              6*CW*(CA*Mf2(4,Gen4))
               dup9 = USf(Sfe4,1,3,Gen4)*
     -              (3*CW*(dup3*Mf(3,Gen4)*USfC(Sfe4,2,3,Gen4)) + 
     -              dup7*USfC(Sfe4,1,3,Gen4)) + 
     -              USf(Sfe4,2,3,Gen4)*
     -              (3*CW*(dup1*Mf(3,Gen4)*USfC(Sfe4,1,3,Gen4)) + 
     -              2*dup5*USfC(Sfe4,2,3,Gen4))
               dup10 = USf(Sfe4,1,4,Gen4)*
     -              (3*CW*(dup4*Mf(4,Gen4)*USfC(Sfe4,2,4,Gen4)) + 
     -              dup8*USfC(Sfe4,1,4,Gen4)) + 
     -              USf(Sfe4,2,4,Gen4)*
     -              (3*CW*(dup2*Mf(4,Gen4)*USfC(Sfe4,1,4,Gen4)) - 
     -              2*dup6*USfC(Sfe4,2,4,Gen4))

               loopa1 = loopa1 + Alfa/(24.D0*CW*MW*Pi)*
     -              ((q12-q22)*(CB*dup9*EL*
     -              (2*Bval(bb0,iint3(Gen4,Sfe4)) + 
     -              4*Bval(bb1,iint3(Gen4,Sfe4))) + 
     -              dup10*EL*SB*
     -              (-Bval(bb0,iint4(Gen4,Sfe4)) - 
     -              2*Bval(bb1,iint4(Gen4,Sfe4)))))/(CB*q32*SB*SW)
            enddo
         enddo

** Now we have contributions from (s)leptons
         if (sector .ne. 1) then
            do Gen4 = 1,3
               loopa1 = loopa1 + Alfa/(2.D0*MW*Pi)*((q12-q22)*CA*EL*
     -              (Bval(bb0,iintL1(Gen4))*Mf2(2,Gen4) + 
     -              2*Bval(bb1,iintL1(Gen4))*Mf2(2,Gen4)))/(CB*q32*SW)
            enddo

            do Sfe4 = 1,2
               do Gen4 = 1,3
                  lup1 = MUEC*SA - CA*Af(2,Gen4,Gen4)
                  lup2 = MUE*SA - CA*AfC(2,Gen4,Gen4)
                  lup3 = MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen4))
                  lup4 = MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                 2*CW*(CA*Mf2(2,Gen4))
                  lup5 = USf(Sfe4,1,2,Gen4)*
     -                 (CW*(lup2*Mf(2,Gen4)*USfC(Sfe4,2,2,Gen4)) + 
     -                 lup4*USfC(Sfe4,1,2,Gen4)) + 
     -                 USf(Sfe4,2,2,Gen4)*
     -                 (CW*(lup1*Mf(2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 
     -                 2*lup3*USfC(Sfe4,2,2,Gen4))

                  loopa1 = loopa1 + 
     -                 Alfa/(8.D0*CW*MW*Pi)*
     -                 ((q12-q22)*lup5*EL*
     -                 (Bval(bb0,iintL2(Gen4,Sfe4)) + 
     -                 2*Bval(bb1,iintL2(Gen4,Sfe4))))/(CB*q32*SW)
               enddo
            enddo
         end if


** Lastly, the contributions from SM-type bosons and their superpartners
         if (sector .ge. 3) then
            loopa1 = loopa1 + (-((Alfa*MW)/(8.D0*CW2*Pi)*((q12-q22)*EL*
     -           Bval(bb0,iintB1)*(-(C2B*CAB) + 2*CW2*CBA))/SW) - (Alfa*
     -           MW)/(4.D0*CW2*Pi)*((q12-q22)*EL*Bval(bb1,iintB1)*
     -           (-(C2B*CAB) + 2*CW2*CBA))/SW - (Alfa*MW)/(8.D0*CW2*Pi)*
     -           ((q12-q22)*EL*Bval(bb0,iintB2)*(C2B*CAB + 6*CW2*CBA))/
     -           SW-(Alfa*MW)/(4.D0*CW2*Pi)*((q12-q22)*EL*
     -           Bval(bb1,iintB2)*(C2B*CAB + 6*CW2*CBA))/SW)/q32

** chargino / neutralino contribution
            if (sector .ge. 5) then
               do Cha4 = 1,2
                  bup1 = CA*(UCha(Cha4,2)*VCha(Cha4,1) + 
     -                 UChaC(Cha4,2)*VChaC(Cha4,1)) + 
     -                 SA*(UCha(Cha4,1)*VCha(Cha4,2) + 
     -                 UChaC(Cha4,1)*VChaC(Cha4,2))
                  
                  loopa1 = loopa1 + 
     -                 Alfa/(2.D0*Pi)*((q12-q22)*bup1*EL*
     -                 (Bval(bb0,iintB3(Cha4))*MCha(Cha4) + 
     -                 2*Bval(bb1,iintB3(Cha4))*MCha(Cha4)))/
     -                 (q32*SW*Sqrt(2.D0))
               enddo
            end if
         end if

** End of calculation of h2 formfactors
      end if


** Lastly, calculation of A0 formfactors
      if (higgno .eq. 3) then

** (S)quark contributions
         do Sfe4 = 1,2
            do Gen4 = g4ini, 3
               dup1 = -((MUEC*TB + Af(3,Gen4,Gen4))*USf(Sfe4,2,3,Gen4)*
     -              USfC(Sfe4,1,3,Gen4)) + 
     -              (MUE*TB + AfC(3,Gen4,Gen4))*USf(Sfe4,1,3,Gen4)*
     -              USfC(Sfe4,2,3,Gen4)
               dup2 = -((MUEC + TB*Af(4,Gen4,Gen4))*USf(Sfe4,2,4,Gen4)*
     -              USfC(Sfe4,1,4,Gen4)) + 
     -              (MUE + TB*AfC(4,Gen4,Gen4))*USf(Sfe4,1,4,Gen4)*
     -              USfC(Sfe4,2,4,Gen4)

               loopa1 = loopa1 + (Alfa*cI)/(8.D0*MW*Pi)*
     -              ((q12-q22)*EL*(dup1*
     -              (-2*Bval(bb0,iint3(Gen4,Sfe4)) - 
     -              4*Bval(bb1,iint3(Gen4,Sfe4)))*Mf(3,Gen4) + 
     -              dup2*TB*(Bval(bb0,iint4(Gen4,Sfe4)) + 
     -              2*Bval(bb1,iint4(Gen4,Sfe4)))*Mf(4,Gen4)))/
     -              (q32*SW*TB)
            enddo
         enddo


** (S)lepton contributions
         if (sector .ne. 1) then
            do Sfe4 = 1,2
               do Gen4 = 1,3
                  lup1 = -((MUEC + TB*Af(2,Gen4,Gen4))*
     -                 USf(Sfe4,2,2,Gen4)*USfC(Sfe4,1,2,Gen4)) + 
     -                 (MUE + TB*AfC(2,Gen4,Gen4))*USf(Sfe4,1,2,Gen4)*
     -                 USfC(Sfe4,2,2,Gen4)

                  loopa1 = loopa1 + (Alfa*cI)/(8.D0*MW*Pi)*((q12-q22)*
     -                 lup1*EL*(Bval(bb0,iintL2(Gen4,Sfe4)) + 2*
     -                 Bval(bb1,iintL2(Gen4,Sfe4)))*Mf(2,Gen4))/(q32*SW)
               enddo
            enddo
         end if


** Contribution from SM-type bosons and their superpartners
         if (sector .ge. 5) then
            do Cha4 = 1,2
               bup1 = SB*(UCha(Cha4,2)*VCha(Cha4,1) - 
     -              UChaC(Cha4,2)*VChaC(Cha4,1)) + 
     -              CB*(UCha(Cha4,1)*VCha(Cha4,2) - 
     -              UChaC(Cha4,1)*VChaC(Cha4,2))
               
               loopa1 = loopa1 + (Alfa*cI)/(2.D0*Pi)*
     -              ((q12-q22)*bup1*EL*
     -              (Bval(bb0,iintB3(Cha4)) + 2*Bval(bb1,iintB3(Cha4)))*
     -              MCha(Cha4))/(q32*SW*Sqrt(2.D0))
            enddo
         end if         

** End of A0 formfactor calculation
      end if


      a1gamh = loopa1


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the formfactors involved in the tree level higgs
*** propagator corrections

      subroutine hgs_prop(a1propW,a1propZ,treeW,treeZ,Hgs,debug,higgno)      

      implicit none


** Higgs propagator-correcting matrix
      double complex Hgs(3,3)


** These are the tree level coupling, the propagator multiplicative factor
** and the propagator contribution itself.  The first entry in each array is
** for the heavy higgs, and the second entry is for the light higgs.
      double complex treeW(3), treeZ(3)
      double complex a1propW, a1propZ


** This parameter defines which of the neutral higgs bosons is being considered.
** It is set in input.dat, line 12
      integer higgno

** Debugging parameter
      integer debug


** Clearing all factors 
      a1propW = (0D0,0D0) 
      a1propZ = (0D0,0D0)


      if (higgno .eq. 1) then
         a1propW = (Hgs(1,1)-1)*treeW(1) + Hgs(1,2)*treeW(2)
         a1propZ = (Hgs(1,1)-1)*treeZ(1) + Hgs(1,2)*treeZ(2)
      else if (higgno .eq. 2) then
         a1propW = Hgs(2,1)*treeW(1) + (Hgs(2,2)-1)*treeW(2)
         a1propZ = Hgs(2,1)*treeZ(1) + (Hgs(2,2)-1)*treeZ(2)
      else if (higgno .eq. 3) then
         a1propW = Hgs(3,1)*treeW(1) + Hgs(3,2)*treeW(2)
         a1propZ = Hgs(3,1)*treeZ(1) + Hgs(3,2)*treeZ(2)
      end if


      if (debug .ne. 0) then
         write(*,*)"  "
         write(*,*)"Higgs propagator corrections ..."
         
         write(*,*)"WW-higgs prop a1 =", a1propW
         write(*,*)"ZZ-higgs prop a1 =", a1propZ
      end if


      end 


********************************************************************************
