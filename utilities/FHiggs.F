********************************************************************************
********************************************************************************
***                                                                          ***
*** FHiggs.F                                                                 ***
*** 18th May 2011                                                            ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This file contain all of the subroutines that talk to directly to        ***
*** FeynHiggs.  It's included in the make iff we compile linked to FeynHiggs.***
***                                                                          ***
*** IMPORTANT NOTE: FeynHiggs changed quite a lot between version 2.6.x and  ***
***                 2.7.x and 2.8.x.  The configure file should work out     ***
***                 which version of FeynHiggs the code is linked to, and    ***
***                 use the appropriate subroutine calls, but future changes ***
***                 in FeynHiggs could cause problems, and older versions    ***
***                 have not been tested.                                    ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** This subroutine deals with all the FeynHiggs parameter settings: it calls
*** the flag setting and Higgs calculation etc. routines.

      subroutine FHpara(t1_CplxApprox,MSusy,scalefactor,DeltaMB,debug)

      implicit none

** Soft-SUSY breaking parameter in the sfermion sector
      double precision MSusy

** Scale information for FeynHiggs
      double precision scalefactor

** Delta MB from FeynHiggs
      double complex DeltaMB

** Note that t1_CplxApprox is set according to whether we're using complex 
** parameters in the MSSM
      integer mssmpart, fieldren, tanbren, p2approx
      integer looplevel, loglevel, t1_running_mt, t1_bot_resum, t1_CplxApprox
      parameter (mssmpart = 4)
      parameter (fieldren = 0)
      parameter (tanbren = 0)
      parameter (p2approx = 0)
      parameter (looplevel = 2)
      parameter (loglevel = 2)
      parameter (t1_running_mt = 1)
      parameter (t1_bot_resum = 1)


      double precision AlfasMZ, invAlfa
      double precision CKMlambda, CKMA, CKMrho, CKMeta
      parameter (AlfasMZ = -1D0)
      
      double precision Qtau, Qt, Qb
      parameter (Qtau = 0)
      parameter (Qt = 0)
      parameter (Qb = 0)

      integer error, debug


#include "SLHA.h"
      double complex slhadata(nslhadata)

#include "mssm.inc"
#include "koppln.inc"
#include "global.inc"


      if (debug .eq. 2) then
         call FHSetDebug(3)
      end if


* setting Flags: this is basically FH initialisation
#if FHVERSION >= 21200
      call FHSetFlags(error,mssmpart,fieldren,tanbren,higgsmix, 
     -     p2approx,looplevel,loglevel,t1_running_mt,t1_bot_resum,t1_CplxApprox)
#else
      call FHSetFlags(error,mssmpart,fieldren,tanbren,higgsmix, 
     -     p2approx,looplevel,t1_running_mt,t1_bot_resum,t1_CplxApprox)
#endif
      if (error .ne. 0) then
         write(*,*)'Problem with FeynHiggs set flags'
         stop
      end if

      write(*,*)'  '
      write(*,*)'FeynHiggs Flags ...'
      write(*,*)'mssmpart =', mssmpart
      write(*,*)'fieldren =', fieldren
      write(*,*)'tanbren =', tanbren
      write(*,*)'higgsmix =', higgsmix
      write(*,*)'p2approx =', p2approx
      write(*,*)'looplevel =', looplevel
      write(*,*)'loglevel =', loglevel
      write(*,*)'t1_running_mt =', t1_running_mt
      write(*,*)'t1_bot_resum =', t1_bot_resum
      write(*,*)'t1_CplxApprox =', t1_CplxApprox


** Now we set the parameters: either via SLHA or from the input values
      if (slha_switch) then

         call SLHARead(error, slhadata, slhaFILE, 1)
         if( error .ne. 0 ) then
            write(*,*)'Problem reading SLHA file'
            stop
         end if

#if FHVERSION >= 20800
         call FHSetSLHA(error, slhadata)
#else
         call FHSetSLHA(error, slhadata, scalefactor)
#endif
         if( error .ne. 0 ) then
            write(*,*)'Problem setting SLHA parameters'
            stop
         end if

      else

         CKMlambda = -1D0
         CKMA = -1D0
         CKMrho = -1D0
         CKMeta = -1D0
         
         invAlfa = 1D0/Alfa

#if FHVERSION >= 20800

         call FHSetSMPara(error,
     &        invAlfa, AlfasMZ, GF,
     &        ME, MsU, MD, MM, MC, MS, ML, MB, MW, MZ,
     &        CKMlambda, CKMA, CKMrho, CKMeta)
#else
         call FHSetSMPara(error,
     &        invAlfa, AlfasMZ, GF,
     &        MS, MC, MB, MW, MZ,
     &        CKMlambda, CKMA, CKMrho, CKMeta)
#endif
         if( error .ne. 0 ) then
            write(*,*)'Problem with FeynHiggs, set SM Para'
            stop
         end if

         call FHSetPara(error, scalefactor,
     &        MT, TB, MA0, MHp,
     &        M3SL, M3SE, M3SQ, M3SU, M3SD,
     &        M2SL, M2SE, M2SQ, M2SU, M2SD,
     &        M1SL, M1SE, M1SQ, M1SU, M1SD,
     &        MUE,
     &        Af(2,3,3), Af(3,3,3), Af(4,3,3),
     &        Af(2,2,2), Af(3,2,2), Af(4,2,2),
     &        Af(2,1,1), Af(3,1,1), Af(4,1,1),
     &        M_1, M_2, M_3,
     &        Qtau, Qt, Qb)
         if( error .ne. 0 ) then
            write(*,*)'Problem with FeynHiggs, set Para'
            stop
         end if
      end if


** Getting the FeynHiggs-calculated masses and mixings of the SUSY particles
      call getPara(MSf,USf,UCha,VCha,ZNeu,MCha,MNeu,DeltaMB)


** Calculating the Higgs masses and mixings
      call HiggsCorrFH(SAeff,MHiggs,Hgs,ZHiggs)


** Calculating the Higgs widths and branching ratios
      call BRs_FeynH(debug)


      end


*******************************************************************************
*******************************************************************************

*** This subroutine retrieves all of the useful mssm parameters calculated by 
*** setPara

      subroutine getPara(MSf,USf,UCha,VCha,ZNeu,MCha,MNeu,DeltaMB)

      implicit none

      integer nmfv
      double precision MSf(2,4,3), MASf(6,5), MCha(2), MNeu(4)
      double complex USf(2,2,4,3), UASf(6,6,5)
      double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
      double complex DeltaMB
      double precision MGl

      double precision MHtree(4), SAtree
#if FHVERSION >= 20905
      double precision MrSf(2,5,3)
      double complex UrSf(2,2,5,3)
#endif

      integer error

** do loop variables used to set sfermion parameters
      integer type, gen, sfe1, sfe2


#if FHVERSION >= 20700

#if FHVERSION >= 20905
** For FeynHiggs versions 2.9.5 and later
      call FHGetPara(error, nmfv, MrSf, UrSf, MASf, UASf, MCha, UCha, VCha, MNeu, 
     &     ZNeu, DeltaMB, MGl, MHtree, SAtree)
#elif FHVERSION >= 20800
** For FeynHiggs versions 2.8.x and 2.9.[0-4]
      call FHGetPara(error, nmfv, MASf, UASf, MCha, UCha, VCha, MNeu, 
     &     ZNeu, DeltaMB, MGl, MHtree, SAtree)
#elif FHVERSION >= 20700
** For FeynHiggs versions 2.7.x
      call FHGetPara(error, nmfv, MASf, UASf, MCha, UCha, VCha, MNeu, 
     &     ZNeu, DeltaMB, MGl)
#endif

** Setting sfermion masses and mixings from the extended squark masses and 
** mixings provided by FeynHiggs
      do type = 1, 4
         do gen = 1, 3
            do sfe1 = 1, 2
               MSf(sfe1,type,gen) = MASf(gen+3*(sfe1-1),type) 
               do sfe2 = 1, 2
                  USf(sfe1,sfe2,type,gen) = 
     -                 UASf(gen+3*(sfe1-1), gen+3*(sfe2-1),type)
               end do
            end do
         end do
      end do


#else
** For FeynHiggs versions 2.6.x, the following code must be used
      call FHGetPara(error, nmfv, MSf, USf, MASf, UASf, MCha, UCha, 
     &     VCha, MNeu, ZNeu, DeltaMB, MGl)
#endif
      
      if( error .ne. 0 ) then
         write(*,*)'Sorry, there is a problem with FH getPara'
         stop
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the higgs sector of the mssm

      subroutine HiggsCorrFH(SAeff,MHiggs,Hgs,ZHiggs)

      implicit none

      integer error
      double precision MHiggs(4), SAeff
      double complex cSAeff, UHiggs(3,3), ZHiggs(3,3)

** When prop = 1, the U matrix is used, when prop = 2, the Z matrix is used.
** The Z matrix is hard-wired
      integer prop
      parameter (prop = 2)

** The do loop variables     
      integer h, j

** The propagator correction matrix: chosen by the parameter 'prop'
      double complex Hgs(3,3)


      call FHHiggsCorr(error,MHiggs,cSAeff,UHiggs,ZHiggs)
      SAeff = dble(cSAeff)

      if (error .ne. 0) then
         write(*,*)'Problem with FeynHiggs calculation of Higgs sector'
         stop
      end if


** Setting Higgs propagator corrections
      if (prop .eq. 1) then
         do h = 1,3
            do j = 1,3
               Hgs(h,j) = UHiggs(h,j)
            end do
         end do 
      else if (prop .eq. 2) then
         do h = 1,3
            do j = 1,3
               Hgs(h,j) = ZHiggs(h,j)
            end do
         end do 
      else
         do h = 1,3
            do j = 1,3
               Hgs(h,j) = (0D0, 0D0)
            end do
         end do
         Hgs(1,1) = (1D0,0D0)
         Hgs(2,2) = (1D0,0D0)
      end if

      
      end


********************************************************************************
********************************************************************************

*** This subroutine uses FeynHiggs to set the branching ratios and widths of the
*** Higgs

*** NOTE: As of 26-04-2011, FeynHiggs had a problem with calculating the
*** decay widths if it was linked dynamically (using the flag -fPIC).  Thus,
*** if VBFNLO is dynamically linked, we don't use FeynHiggs to calculate the
*** branching ratios, but instead either take values from a SLHA file or from
*** an internal calculation (in utilities/koppln.F).

      subroutine BRs_FeynH(debug)

      implicit none
   
*sophy: added FeynHiggs - Higgs widths
#include "FHCouplings.h"
#include "global.inc"

      integer i

** FeynHiggs outputs:
      integer error
      double complex couplings(ncouplings)
      double complex couplingsms(ncouplingsms)
      double precision gammas(ngammas)
      double precision gammasms(ngammasms)
      integer excl 

** Total widths and branching ratios: VBFNLO
      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT


      integer debug


** Function that determines whether a parameter needs to be reset here.
      logical replace
      external replace


**clearing values of gammas
      gammas = 0D0
      gammasms = 0D0


** Using FH to calculate the widths and branching ratios of the Higgs
** Only SM-type decays are included for now.  See notebook 10-1-11 for indices

** For versions of FeynHiggs 2.6.x
#if FHVERSION >= 20600 && FHVERSION < 20700
      call FHCouplings(error,couplings,couplingsms,
     -     gammas,gammasms,excl)

      if (model .eq. 2) then
** susy ..
         if (replace(BHZZ,0)) BHZZ = BR(H0VV(HiggsType,3))
         if (replace(BHWW,0)) BHWW = BR(H0VV(HiggsType,4))
         if (replace(BHBB,0)) BHBB = BR(H0FF(HiggsType,4,3))
         if (replace(BHCC,0)) BHCC = BR(H0FF(HiggsType,3,2))
         if (replace(BHTT,0)) BHTT = BR(H0FF(HiggsType,3,3))
         if (replace(BHtau,0)) BHtau = BR(H0FF(HiggsType,2,3))
         if (replace(BHmu,0)) BHmu= BR(H0FF(HiggsType,2,2))
         if (replace(BHgam,0)) BHgam = BR(H0VV(HiggsType,1))
         if (replace(BHgamZ,0)) BHgamZ = BR(H0VV(HiggsType,2))

         if (replace(XGH,0)) XGH = GammaTot(HiggsType)
      else
** sm ..
         if (replace(BHZZ,0)) BHZZ = BRSM(H0VV(HiggsType,3))
         if (replace(BHWW,0)) BHWW = BRSM(H0VV(HiggsType,4))
         if (replace(BHBB,0)) BHBB = BRSM(H0FF(HiggsType,4,3))
         if (replace(BHCC,0)) BHCC = BRSM(H0FF(HiggsType,3,2))
         if (replace(BHTT,0)) BHTT = BRSM(H0FF(HiggsType,3,3))
         if (replace(BHtau,0)) BHtau = BRSM(H0FF(HiggsType,2,3))
         if (replace(BHmu,0)) BHmu= BRSM(H0FF(HiggsType,2,2))
         if (replace(BHgam,0)) BHgam = BRSM(H0VV(HiggsType,1))
         if (replace(BHgamZ,0)) BHgamZ = BRSM(H0VV(HiggsType,2))

         if (replace(XGH,0)) XGH = GammaSMTot(HiggsType)
      end if

#else

      call FHCouplings(error,couplings,couplingsms,
     -     gammas,gammasms,1)
** FeynHiggs 2.7.x or 2.8.x

      if (model .eq. 2) then
** susy ..         
         if (replace(BHZZ,0)) BHZZ = BR(H0VV(HiggsType,3))
         if (replace(BHWW,0)) BHWW = BR(H0VV(HiggsType,4))
         if (replace(BHBB,0)) BHBB = BR(H0FF(HiggsType,4,3,3))
         if (replace(BHCC,0)) BHCC = BR(H0FF(HiggsType,3,2,2))
         if (replace(BHTT,0)) BHTT = BR(H0FF(HiggsType,3,3,3))
         if (replace(BHtau,0)) BHtau = BR(H0FF(HiggsType,2,3,3))
         if (replace(BHmu,0)) BHmu= BR(H0FF(HiggsType,2,2,2))
         if (replace(BHgam,0)) BHgam = BR(H0VV(HiggsType,1))
         if (replace(BHgamZ,0)) BHgamZ = BR(H0VV(HiggsType,2))
         if (replace(BHgg,0)) BHgg = BR(H0VV(HiggsType,5))
         
         if (replace(XGH,0)) XGH = GammaTot(HiggsType)

         write(*,*)'BHZZ =', BHZZ
         write(*,*)'BHWW =', BHWW
         write(*,*)'BHBB =', BHBB
         write(*,*)'BHCC =', BHCC
         write(*,*)'BHTT =', BHTT
         write(*,*)'BHtau =', BHtau
         write(*,*)'BHmu =', BHmu
         write(*,*)'BHgam =', BHgam
         write(*,*)'BHgamZ =', BHgamZ
         write(*,*)'BHgg =', BHgg

      else
** sm ..         
         if (replace(BHZZ,0)) BHZZ = BRSM(H0VV(HiggsType,3))
         if (replace(BHWW,0)) BHWW = BRSM(H0VV(HiggsType,4))
         if (replace(BHBB,0)) BHBB = BRSM(H0FF(HiggsType,4,3,3))
         if (replace(BHCC,0)) BHCC = BRSM(H0FF(HiggsType,3,2,2))
         if (replace(BHTT,0)) BHTT = BRSM(H0FF(HiggsType,3,3,3))
         if (replace(BHtau,0)) BHtau = BRSM(H0FF(HiggsType,2,3,3))
         if (replace(BHmu,0)) BHmu= BRSM(H0FF(HiggsType,2,2,2))
         if (replace(BHgam,0)) BHgam = BRSM(H0VV(HiggsType,1))
         if (replace(BHgamZ,0)) BHgamZ = BRSM(H0VV(HiggsType,2))
         if (replace(BHgg,0)) BHgg = BRSM(H0VV(HiggsType,5))
         
         if (replace(XGH,0)) XGH = GammaSMTot(HiggsType)
      end if
#endif

      if( error .ne. 0 ) then
         write(*,*)'Sorry, there was a problem with FHCouplings'
         stop
      end if


      end


********************************************************************************
