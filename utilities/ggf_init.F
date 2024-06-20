* higgs_mix= 2  ! Mixing between MSSM Higgs bosons                                                                      !                     *               !  according to ZHiggs     
* higgs_mix=1   ! Mixing cp-odd and -even, user defined
* higgs_mix= 0  ! No higgs mixing: default case   
* model=1       ! SM Higgs
* model=2       ! MSSM, no higgs mixing       
* HiggsType= 0,1,2,3  0:SM      
***********************************************************************
      SUBROUTINE init_gluonfusion
         use globalvars, only: lglobalprint
          use readinput
*----------------------------------------------------------------------
c... author   : Michael Kubocz
c... modified : 28.08.2013
c... author   : Francisco Campanario 
c... modified : 12.01.2017
***********************************************************************
      use globalvars, only: lglobalprint
      IMPLICIT NONE

      logical MGdebug,lgauge2
      common/MGdebugC/MGdebug,lgauge2
      logical MGdebug1,MGPen,MGHex,PRINT_AMP
      common/MGdebugC1/MGdebug1,MGPen,MGHex,PRINT_AMP
      real*8 series_t,series_p, epsmt
      complex*16 ffactorH, ffactorA
c     series_t and series_p: To be found in MADGRAPH in SOURCE/MODEL/couplings.f
c     Define with the same names.   
      integer icall1
      data icall1/0/
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/mssm.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/amplitudes/ggf/ggf_process.inc"
#include "VBFNLO/amplitudes/ggf/ggf_default.inc"

      CALL loadfile("ggflo.dat",.true.)
c... Loop particles choice
      CALL read_int("Q_LOOP",Q_loop,3)
c... subprocess choice
      CALL read_Logical("SUBPRQQ",subprQQ,.true.,.true.)
      CALL read_Logical("SUBPRQG",subprQG,.true.,.true.)
      CALL read_Logical("SUBPRGG",subprGG,.true.,.true.)
      CALL read_int("SUBPRQG_PART",sub_QG_part,0,.true.)
      CALL read_int("SUBPRQQ_PART",sub_QQ_part,0,.true.)
c... Mixing between different Higgs bosons
      CALL read_int("HIGGS_MIX",higgs_mix,0)
cfc      if (higgs_mix.ne.2) then
      if (higgs_mix.eq.1) then
         CALL read_real("CP_EVEN_MOD",cp_even_mod,1d0)
         CALL read_real("CP_ODD_MOD",cp_odd_mod,1d0)

         select CASE (procID)
         CASE(GFHjj, GFHjj_WW, GFHjj_ZZ_ll,GFHjj_ZZ_lnu,GFHjjj)
         CASE DEFAULT
            write(*,*)'  '
            write(*,*)'Sorry, we can not yet include higgs decays'
            write(*,*)'into fermions or photons from gluon fusion with '
            write(*,*)'user-defined mixing between CP-odd and '
            write(*,*)'CP-even Higgs bosons'
            stop
         END select
      end if
      call read_logical("EFF_TH_CORR",eff_th_corr,.true.,.true.) 
      call read_logical("EFF_TH_CORR_TAYLOR",eff_th_corr_taylor,.false.,.true.) 
      call read_logical("DBLE_SWITCH",dble1,.true.,.true.)  
      call read_logical("QUAD_SWITCH",QUAD,.false.,.true.)  
      call read_Real("GAUGELIMIT",gaugelimit,0.01d0,.true.)
      call read_logical("GAUGE_SWITCH",lgauge2,.false.,.true.)  
      call read_logical("PENTAGONSTO0",MGPen,.false.,.true.)  
      call read_logical("HEXAGONSTO0",MGHex,.false.,.true.) 
      
      call read_logical("HEX_EVALUATE",hexeva,.true.,.true.)  
      call read_logical("PEN_EVALUATE",peneva,.true.,.true.) 
      call read_logical("BOX_EVALUATE",boxeva,.true.,.true.) 

      call read_logical("MG_PRINT",MGdebug1,.false.,.true.)  
      call read_logical("MG_DEBUG",MGdebug,.false.,.true.)        
      call read_logical("MG_UU",mguu,.false.,.true.) 
      call read_logical("MG_UUB",mguub,.false.,.true.) 
      call read_logical("MG_GU",mggu,.false.,.true.) 
      call read_logical("MG_EQUAL_FLAVOUR",equal_f,.false.,.true.) 

      call read_logical("PRINT_AMP",PRINT_AMP,.false.,.true.) 
      call read_Real("GRAMCTO0",GRAMCTO0,1d-6,.true.)
      call read_Real("GRAMDTO0",GRAMDTO0,1d-6,.true.)

      call printnfl(.false.) ! If true indicates that a W is attached to quark line. This is relevant of VBF processes.
 

      CALL closefile
 
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------
*----------------------------------------------------------------------
      tan_Beta  = TB 
      sin_Beta  = SB
      cos_Beta  = CB
      sin_Alpha = SAtree
      cos_Alpha = CAtree
!**********************************************************************
      if (lglobalprint)  write(6,*)'--------- PROCESS PARAMETERS IN GLUON FUSION ---------'
      if (lglobalprint)  write(6,*)''
      if (abs(xmh-2d0*m_top) .le. 1d0) then
         if (lglobalprint) then
        write(6,*)' -->                  ATTENTION                  <--'
        write(6,*)' -->     threshold effects at m_H = 2*m_top      <--'
        write(6,*)' -->    please use a Higgs mass that fulfills    <--'
        write(6,*)' -->          | m_H - 2*m_top | > 1 GeV          <--'
        write(6,*)' -->                  ATTENTION                  <--'
        write(6,*)''
        endif
      endif
      if (lglobalprint) write(6,*)'------- FURTHER PROCESS INFO  -------'
      if (lglobalprint) then
         if(Q_loop.eq.1) write(6,"(T4,A,T22,A,F16.4)")"m_top","= ",m_top
         if(Q_loop.eq.2) then 
            write(6,"(T4,A,T22,A,F9.4)")"m_bottom(m_Higgs)","= ",BMASS_H
            write(6,"(T4,A,T22,A,F9.4)")"m_bottom pole mass in loops"
     &         ,"=" ,m_bottom
         endif
         if(Q_loop.eq.3) then 
            write(6,"(T4,A,T22,A,F16.4)")"m_top","= ",m_top
            write(6,"(T4,A,T22,A,F9.4)")"m_bottom(m_Higgs)","= ",BMASS_H
            write(6,"(T4,A,T22,A,F9.4)")"m_bottom pole mass in loops"
     &         ,"=" ,m_bottom
         endif
         
         print*,''
         write(6,*)'Gram det in Cij ->', GRAMCTO0
         write(6,*)'Gram det in Dij ->', GRAMDTO0
         
         write(6,*)'Gauge limit ->',  Gaugelimit
         write(6,*)'Evaluate Dble ->', dble1
         write(6,*)'Evaluate QUAD ->', QUAD
         
         print*,''
         
         write(6,*) '------------------- PROCESS(ES) ----------------'
         
         if (subprQQ) write(6,*) 'QUARK - QUARK SCATTERING '
         if (subprQG) write(6,*) 'QUARK - GLUON SCATTERING ' 
         if (subprGG) write(6,*) 'GLUON - GLUON SCATTERING '
         if (sub_QG_PART.ne.0) THEN
            IF (sub_QG_PART.EQ.1) THEN
               WRITE(6,*) 'type: gg :> q~q g' 
            elseif (sub_QG_PART.EQ.2) THEN
               WRITE(6,*) 'type: gq :> gg q ' 
            elseif (sub_QG_PART.EQ.3) THEN
               WRITE(6,*) 'type: qq :> ggg '
            else
               WRITE(6,*) 'wrong type:'
               stop
            endif
         endif
         

         if (sub_QQ_PART.ne.0) THEN
            IF (sub_QQ_PART.EQ.1) THEN
               WRITE(6,*) 'type: qq :> q~q g' 
            elseif (sub_QQ_PART.EQ.2) THEN
               WRITE(6,*) 'type: gq :> q~q q~q ' 
            else
               WRITE(6,*) 'wrong type:'
               stop
            endif
         endif
      endif

***********************************************************************

      if (higgs_mix .eq. 0) then   ! No higgs mixing: default case
         if (model.eq.1) then      ! SM Higgs
            if(Q_loop.eq.0) then
               if (lglobalprint) then
                  print*,''
                  write(6,*)'cp-even SM Higgs in the mtop->infinity limit'
               endif
               Heff          = .true.
cccccccccccccccccccccccccccccccccccccccccccccccccc
cMK   Form factor correction to the eff-th. vertex (see hep-ph/0503172, p.89)
               ffactorH = dcmplx(0d0,0d0) 
               if(eff_th_corr) then
                  epsmt = xmh**2/(4d0*m_top**2)
                  if(xmh.lt.2d0*m_top) then
                     if(eff_th_corr_taylor) then
cFC   To agree with Madgraph 4, add terms of the expansion
                        series_t = 1d0 + epsmt*7d0/30d0 + epsmt**2*2d0/21d0 
     &                       + epsmt**3*26d0/525d0
                        ffactorH = dcmplx(series_t,0d0)
                        if (lglobalprint)  Print*,"CP-even coupling correction from 1/M_TOP series expansion of FF",series_t
                     else
                        series_t = 2d0*(epsmt+(epsmt-1d0)
     &                       *(DASIN(DSQRT(epsmt)))**2)/(epsmt**2)
                        ffactorH = dcmplx(series_t,0d0)*3d0/4d0
                        if (lglobalprint)  Print*,"CP-even coupling correction with full FF",ffactorH
                     endif
                  elseif(xmh.gt.2d0*m_top) then
                     print*,'M_HIGGS > 2*M_TOP'
                     if(eff_th_corr_taylor) then
                        if (lglobalprint) print*,'1/M_TOP Series switched off'
                        ffactorH = dcmplx(1d0,0d0)
                     else
                        ffactorH =2d0*(epsmt+(epsmt-1d0)*(-0.25d0*
     &                       (DLOG((1d0+DSQRT(1d0-1d0/epsmt))/(1d0
     &                       -DSQRT(1d0-1d0/epsmt)))-(0d0,1d0)
     &                       *3.141592653589793d0)**2))/(epsmt**2)
     &                       *3d0/4d0
                        if (lglobalprint)  Print*,"CP-even coupling correction with full FF",ffactorH
                     endif
                  endif
               else
                  ffactorH = dcmplx(1d0,0d0)
                  if (lglobalprint)  Print*,"CP-even coupling with no FF corrections",ffactorH
               endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               cplModHeff    = 1d0/3d0*ffactorH
            else if(Q_loop.eq.1) then
               if (lglobalprint) then
                  print*,''
                  write(6,*) 'cp-even SM Higgs (coupling via top-loop) '
               endif
               HloopTop      = .true.
               cplModHtop    = m_top**2
            else if(Q_loop.eq.2) then
               if (lglobalprint) then
                  print*,''
                  write(6,*) 'cp-even SM Higgs (coupling via bottom-loop) '
               endif

               HloopBottom   = .true.
               cplModHbottom = BMASS_H**2
            else if(Q_loop.eq.3) then
               if (lglobalprint) then
                  print*,''
                  write(6,*) 'cp-even SM Higgs (coupling via top-loop) '
                  write(6,*) 'cp-even SM Higgs (coupling via bottom-loop) '
               endif
               HloopTop      = .true.
               HloopBottom   = .true.
               cplModHtop    = m_top**2
               cplModHbottom = BMASS_H**2
            endif
         endif

         if(model.eq.2) then    ! MSSM, no higgs mixing
*       if(HiggsType.eq.0) then
*         write(*,*)'When working in the MSSM, please choose'
*         write(*,*)'HTYPE = 1,2,3'
*         STOP
*       endif
            if(HiggsType.eq.3) then  ! CP odd higgs in MSSM, no mixing
               if(Q_loop.eq.0) then

                  if (lglobalprint) then
                     print*,' '
                     write(6,*)'cp-odd MSSM Higgs A'
                     write(6,*)'in the mtop->infinity limit'
                  endif
                  Aeff          = .true.
cccccccccccccccccccccccccccccccccccccccccccccccccc
cMK   Form factor correction to the eff-th. vertex (see hep-ph/0503173, p.92)
               ffactorA = dcmplx(0d0,0d0) 
               if(eff_th_corr) then
                  epsmt = xmh**2/(4d0*m_top**2)
                  if(xmh.lt.2d0*m_top) then
                     if(eff_th_corr_taylor) then
cFC   To agree with Magraph 4, add terms of the expansion
                     series_p = 1d0 + epsmt/3d0  
     &                       + epsmt**2*8d0/45d0 + epsmt**3*4d0/35d0
                     ffactorA = dcmplx(series_p,0d0)
                     if (lglobalprint)  Print*,"CP-odd coupling correction from 1/M_TOP series expansion of FF",series_p
                  else
                     series_p = (2d0*(DASIN(DSQRT(epsmt)))**2)/epsmt
                     ffactorA = dcmplx(series_p,0d0)/2d0
                     if (lglobalprint)  Print*,"CP-odd coupling correction with full FF",ffactorA
                  endif
               elseif(xmh.gt.2d0*m_top) then
                  print*,'M_HIGGS > 2*M_TOP'
                  if(eff_th_corr_taylor) then
                     if (lglobalprint) print*,'1/M_TOP Series switched off'
                     ffactorA = dcmplx(1d0,0d0)
                  else
                     ffactorA =2d0*(-0.25d0*(DLOG((1d0+DSQRT(1d0-1d0
     &                    /epsmt))/(1d0-DSQRT(1d0-1d0/epsmt)))-(0d0,1d0)
     &                    *3.141592653589793d0)**2)/epsmt/2d0
                     if (lglobalprint)  Print*,"CP-odd coupling correction with full FF",ffactorA
                  endif
               endif
            else
               ffactorA = dcmplx(1d0,0d0)
               if (lglobalprint)  Print*,"CP-odd coupling with no FF corrections",ffactorA
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
            cplModAeff    = 0.5d0*ffactorA
                  if (lglobalprint)  print*,'cplModAeff',cplModAeff
               elseif(Q_loop.eq.1) then
                  if (lglobalprint) then
                     print*,' '
                     write(6,*) 'cp-odd MSSM Higgs A:'
                     write(6,*) '(coupling via top-loop) '
                     write(6,*) 'Tan_BETA ', tan_Beta
                  endif
                  AloopTop      = .true.
                  cplModAtop    = (m_top**2)/tan_Beta
               elseif(Q_loop.eq.2) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-odd MSSM Higgs A:'
                     write(6,*) '(coupling via bottom-loop) '
                     write(6,*) 'Tan_BETA ', tan_Beta
                  endif
                  AloopBottom   = .true.
                  cplModAbottom = (BMASS_H**2)*tan_Beta
               elseif(Q_loop.eq.3) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-odd MSSM Higgs A:'
                     write(6,*) '(coupling via top-loop) '
                     write(6,*) '(coupling via bottom-loop) '
                     write(6,*) 'Tan_BETA ', tan_Beta
                  endif
                  AloopTop      = .true.
                  AloopBottom   = .true.
                  cplModAtop    = (m_top**2)/tan_Beta
                  cplModAbottom = (BMASS_H**2)*tan_Beta
               endif
            endif

            if(HiggsType.eq.1) then  ! MSSM, light cp-even h0, no higgs mixing
               if(Q_loop.eq.0) then
                  if (lglobalprint) write(*,*)'process is not available'
                  STOP
               endif
               if(Q_loop.eq.1) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-even MSSM Higgs h:'
                     write(6,*) '(coupling via top-loop) '
                  endif
                  sHloopTop     = .true.
                  cplModsHtop   =  -(m_top**2)*cos_alpha/sin_Beta
               elseif(Q_loop.eq.2) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-even MSSM Higgs h:'
                     write(6,*) '(coupling via bottom-loop) '
                  endif
                  sHloopBottom  = .true.
                  cplModsHbottom= (BMASS_H**2)*sin_alpha/cos_Beta
               elseif(Q_loop.eq.3) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-even MSSM Higgs h:'
                     write(6,*) '(coupling via top-loop) '
                     write(6,*) '(coupling via bottom-loop) '
                  endif
                  sHloopTop     = .true.
                  sHloopBottom  = .true.
                  cplModsHtop   = -(m_top**2)*cos_alpha/sin_Beta
                  cplModsHbottom= (BMASS_H**2)*sin_alpha/cos_Beta 
               endif
            endif

            if(HiggsType.eq.2) then  !  MSSM, cp-even H, no higgs mixing
               if(Q_loop.eq.0) then
                  if (lglobalprint) write(*,*)'process is not available'
                  STOP
               endif
               if(Q_loop.eq.1) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-even MSSM Higgs H:'
                     write(6,*) '(coupling via top-loop) '
                  endif
                  bHloopTop     = .true.
                  cplModbHtop   = -(m_top**2)*sin_alpha/sin_Beta 
               elseif(Q_loop.eq.2) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-even MSSM Higgs H:'
                     write(6,*) '(coupling via bottom-loop) '
                  endif
                  bHloopBottom  = .true.
                  cplModbHbottom= -(BMASS_H**2)*cos_alpha/cos_Beta 
               elseif(Q_loop.eq.3) then
                  if (lglobalprint) then
                     print*,''
                     write(6,*) 'cp-even MSSM Higgs H:'
                     write(6,*) '(coupling via top-loop) '
                     write(6,*) '(coupling via bottom-loop) '
                  endif
                  bHloopTop     = .true.
                  bHloopBottom  = .true.
                  cplModbHtop   = -(m_top**2)*sin_alpha/sin_Beta
                  cplModbHbottom= -(BMASS_H**2)*cos_alpha/cos_Beta 
               endif
            endif
         endif
      

      else if (higgs_mix .eq. 1) then ! Mixing cp-odd and -even, user defined
         if (lglobalprint) then
            write(6,*)'CP-violating model:'
            write(6,*)'Higgs coupling->(1+i*gamma_5)'
            write(*,"(T4,A,T22,A,F9.4)")"cp_even_mod","= ",cp_even_mod
            write(*,"(T4,A,T22,A,F9.4)")"cp_odd_mod","= ",cp_odd_mod
            print*,''
         endif

         if ((procID .ne.4100).or.(procID .ne.4200)) then
            if (lglobalprint) then
               write(*,*)'  '
               write(*,*)'Sorry, we can not yet include higgs decays'
               write(*,*)'from gluon fusion with user-defined mixing'
               write(*,*)'between CP-odd and CP-even Higgs bosons'
            endif
            stop
         end if
         
         if(model.eq.1) then    ! SM: higgs = cp_even*H_SM + cp_odd*A_SM
!                                 note: A_SM is defined as A_MSSM with TB=1
            if(Q_loop.eq.0) then

                if (lglobalprint) write(6,*) 'cp-even Higgs in the mtop->infinity limit'
               Heff          = .true.
cccccccccccccccccccccccccccccccccccccccccccccccccc
cMK   Form factor correction to the eff-th. vertex (see hep-ph/0503172, p.89)
               ffactorH = dcmplx(0d0,0d0) 
               if(eff_th_corr) then
                  epsmt = xmh**2/(4d0*m_top**2)
                  if(xmh.lt.2d0*m_top) then
                     if(eff_th_corr_taylor) then
cFC   To agree with Madgraph 4, add terms of the expansion
                        series_t = 1d0 + epsmt*7d0/30d0 + epsmt**2*2d0/21d0 
     &                       + epsmt**3*26d0/525d0
                        ffactorH = dcmplx(series_t,0d0)
                        if (lglobalprint)  Print*,"CP-even coupling correction from 1/M_TOP series expansion of FF",series_t
                     else
                        series_t = 2d0*(epsmt+(epsmt-1d0)
     &                       *(DASIN(DSQRT(epsmt)))**2)/(epsmt**2)
                        ffactorH = dcmplx(series_t,0d0)*3d0/4d0
                        if (lglobalprint)  Print*,"CP-even coupling correction with full FF",ffactorH
                     endif
                  elseif(xmh.gt.2d0*m_top) then
                     print*,'M_HIGGS > 2*M_TOP'
                     if(eff_th_corr_taylor) then
                        if (lglobalprint) print*,'1/M_TOP Series switched off'
                        ffactorH = dcmplx(1d0,0d0)
                     else
                        ffactorH =2d0*(epsmt+(epsmt-1d0)*(-0.25d0*
     &                       (DLOG((1d0+DSQRT(1d0-1d0/epsmt))/(1d0
     &                       -DSQRT(1d0-1d0/epsmt)))-(0d0,1d0)
     &                       *3.141592653589793d0)**2))/(epsmt**2)
     &                       *3d0/4d0
                        if (lglobalprint)  Print*,"CP-even coupling correction with full FF",ffactorH
                     endif
                  endif
               else
                  ffactorH = dcmplx(1d0,0d0)
                  if (lglobalprint)  Print*,"CP-even coupling with no FF corrections",ffactorH
               endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               cplModHeff    = cp_even_mod*1d0/3d0*ffactorH
                if (lglobalprint) write(6,*) 'cp-odd Higgs in the mtop->infinity limit'
                Aeff          = .true.
cccccccccccccccccccccccccccccccccccccccccccccccccc
cMK   Form factor correction to the eff-th. vertex (see hep-ph/0503173, p.92)
               ffactorA = dcmplx(0d0,0d0) 
               if(eff_th_corr) then
                  epsmt = xmh**2/(4d0*m_top**2)
                  if(xmh.lt.2d0*m_top) then
                     if(eff_th_corr_taylor) then
cFC   To agree with Magraph 4, add terms of the expansion
                     series_p = 1d0 + epsmt/3d0  
     &                       + epsmt**2*8d0/45d0 + epsmt**3*4d0/35d0
                     ffactorA = dcmplx(series_p,0d0)
                     if (lglobalprint)  Print*,"CP-odd coupling correction from 1/M_TOP series expansion of FF",series_p
                  else
                     series_p = (2d0*(DASIN(DSQRT(epsmt)))**2)/epsmt
                     ffactorA = dcmplx(series_p,0d0)/2d0
                     if (lglobalprint)  Print*,"CP-odd coupling correction with full FF",ffactorA
                  endif
               elseif(xmh.gt.2d0*m_top) then
                  print*,'M_HIGGS > 2*M_TOP'
                  if(eff_th_corr_taylor) then
                     if (lglobalprint) print*,'1/M_TOP Series switched off'
                     ffactorA = dcmplx(1d0,0d0)
                  else
                     ffactorA =2d0*(-0.25d0*(DLOG((1d0+DSQRT(1d0-1d0
     &                    /epsmt))/(1d0-DSQRT(1d0-1d0/epsmt)))-(0d0,1d0)
     &                    *3.141592653589793d0)**2)/epsmt/2d0
                     if (lglobalprint)  Print*,"CP-odd coupling correction with full FF",ffactorA
                  endif
               endif
            else
               ffactorA = dcmplx(1d0,0d0)
               if (lglobalprint)  Print*,"CP-odd coupling with no FF corrections",ffactorA
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
               cplModAeff    = cp_odd_mod*0.5d0*ffactorA
            else if(Q_loop.eq.1) then
               if (lglobalprint) write(6,*) 'cp-even Higgs (coupling via top-loop) '
               HloopTop      = .true.
               cplModHtop    = cp_even_mod*m_top**2
               if (lglobalprint) write(6,*) 'cp-odd Higgs (coupling via top-loop) '
               AloopTop      = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)
            else if(Q_loop.eq.2) then
               if (lglobalprint) write(6,*) 'cp-even Higgs (coupling via bottom-loop) '
               HloopBottom   = .true.
               cplModHbottom = cp_even_mod*BMASS_H**2
               if (lglobalprint) write(6,*) 'cp-odd Higgs (coupling via bottom-loop) '
               AloopBottom   = .true.
               cplModAbottom = cp_odd_mod*(BMASS_H**2)
            else if(Q_loop.eq.3) then
               if (lglobalprint) then
                  write(6,*) 'cp-even SM Higgs: (coupling via top-loop)'
                  write(6,*) 'and (coupling via bottom-loop) '
               endif
               HloopTop      = .true.
               HloopBottom   = .true.
               cplModHtop    = cp_even_mod*(m_top**2)
               cplModHbottom = cp_even_mod*(BMASS_H**2)
               if (lglobalprint) then
                  write(6,*) 'cp-odd Higgs: (coupling via top-loop) '
                  write(6,*) 'and (coupling via bottom-loop) '
               endif
               AloopTop      = .true.
               AloopBottom   = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)
               cplModAbottom = cp_odd_mod*(BMASS_H**2)
            endif               ! end of SM Q_loop choice
            
         else !  MSSM: user-defined mixing, higgs = cp_even*h0 + cp_odd*A0

            if(Q_loop.eq.1) then
               if (lglobalprint) write(6,*) 'cp-odd MSSM Higgs A (via top-loop) '
               AloopTop = .true.
               cplModAtop = cp_odd_mod*(m_top**2)/tan_Beta
               if (lglobalprint) write(6,*) 'cp-even MSSM Higgs h (via top-loop) '
               sHloopTop = .true.
               cplModsHtop = -cp_even_mod*(m_top**2)*cos_alpha/
     -              sin_Beta
            else if(Q_loop.eq.2) then
               if (lglobalprint) write(6,*) 'cp-odd MSSM Higgs A (via bottom-loop) '
               AloopBottom   = .true.
               cplModAbottom = cp_odd_mod*(BMASS_H**2)*tan_Beta
               if (lglobalprint) write(6,*) 'cp-even MSSM Higgs h (via bottom-loop) '
               sHloopBottom   = .true.
               cplModsHbottom = cp_even_mod*(BMASS_H**2)*sin_alpha/
     -              cos_Beta
            else if(Q_loop.eq.3) then
               if (lglobalprint) then
                  write(6,*) 'cp-odd MSSM Higgs A (via top-loop) '
                  write(6,*) 'cp-odd MSSM Higgs A (via bottom-loop) '
               endif
               AloopTop      = .true.
               AloopBottom   = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)/tan_Beta
               cplModAbottom = cp_odd_mod*(BMASS_H**2)*tan_Beta
               if (lglobalprint) then
                  write(6,*) 'cp-even MSSM Higgs h (via top-loop) '
                  write(6,*) 'cp-even MSSM Higgs h (via bottom-loop) '
               endif
               sHloopTop      = .true.
               sHloopBottom   = .true.
               cplModsHtop    = -cp_even_mod*(m_top**2)*cos_alpha/
     -              sin_Beta
               cplModsHbottom = cp_even_mod*(BMASS_H**2)*sin_alpha/
     -              cos_Beta
            end if              ! end of Q_loop choice for MSSM
         endif                  ! end of SM/MSSM for user-defined mixing


      else if (higgs_mix .eq. 2) then !  Mixing between MSSM Higgs bosons 
!                                        according to ZHiggs

         if (model .eq. 1) then ! SM
            if (lglobalprint) then
               write(*,*)'Sorry, we can not mix according to the'
               write(*,*)'Z matrix while working in the SM'
            endif
            stop
         else                   ! MSSM
            if (lglobalprint) write(*,*)'Higgs mixing according to Z propagator matrix:'

            if (Q_loop .eq. 0) then
               if (lglobalprint) then
                  write(*,*)'Sorry, we can not mix Higgs bosons '
                  write(*,*)'according to the Z matric with Q_loop = 0'
               endif
               stop
            else if (Q_loop .eq. 1) then ! top loop
               sHloopTop = .true.
               cplModsHtop =  -DBLE(ZHiggs(HiggsType,1))*(m_top**2)*
     -              cos_alpha/sin_Beta
               bHloopTop = .true.
               cplModbHtop = -DBLE(ZHiggs(HiggsType,2))*(m_top**2)*
     -              sin_alpha/sin_Beta 
               AloopTop = .true.
               cplModAtop = DBLE(ZHiggs(HiggsType,3))*(m_top**2)/
     -              tan_Beta
            else if(Q_loop.eq.2) then
               sHloopBottom  = .true.
               cplModsHbottom= DBLE(ZHiggs(HiggsType,1))*(BMASS_H**2)*
     -              sin_alpha/cos_Beta
               bHloopBottom  = .true.
               cplModbHbottom= -DBLE(ZHiggs(HiggsType,2))*(BMASS_H**2)*
     -              cos_alpha/cos_Beta 
               AloopBottom   = .true.
               cplModAbottom = DBLE(ZHiggs(HiggsType,3))*(BMASS_H**2)*
     -              tan_Beta
            else if(Q_loop.eq.3) then
               sHloopTop     = .true.
               cplModsHtop   = -DBLE(ZHiggs(HiggsType,1))*(m_top**2)*
     -              cos_alpha/sin_Beta
               sHloopBottom  = .true.
               cplModsHbottom= DBLE(ZHiggs(HiggsType,1))*(BMASS_H**2)*
     -              sin_alpha/cos_Beta 
               bHloopTop     = .true.
               cplModbHtop   = -DBLE(ZHiggs(HiggsType,2))*(m_top**2)*
     -              sin_alpha/sin_Beta
               bHloopBottom  = .true.
               cplModbHbottom= -DBLE(ZHiggs(HiggsType,2))*(BMASS_H**2)*
     -              cos_alpha/cos_Beta 
               AloopTop      = .true.
               cplModAtop    = DBLE(ZHiggs(HiggsType,3))*(m_top**2)/
     -              tan_Beta
               AloopBottom   = .true.
               cplModAbottom = DBLE(ZHiggs(HiggsType,3))*(BMASS_H**2)*
     -              tan_Beta
            end if              ! end of Q_loop choice
            
         end if                 ! end of SM/MSSM loop
         
c... Mixing between A & H with tan_Beta dependend Yukawa-couplings to Top-
c... and Bottom-quarks (+user defined factors) within an MSSM-scenario
      elseif(higgs_mix.eq.3) then 
         if (lglobalprint) then
         write(6,*)'CP-violating model:'
         write(6,*)'Higgs coupling->(1+i*gamma_5)'
         write(6,*)'Y_top/tan_Beta and Y_bottom*tan_Beta'
         write(*,"(T4,A,T22,A,F9.4)")"tan_Beta","= ",tan_Beta
         write(*,"(T4,A,T22,A,F9.4)")"cot_Beta","= ",1d0/tan_Beta
         write(*,"(T4,A,T22,A,F9.4)")"cp_even_mod","= ",cp_even_mod
         write(*,"(T4,A,T22,A,F9.4)")"cp_odd_mod","= ",cp_odd_mod
         print*,''
         endif
         if(model.eq.1) then  
            if (lglobalprint) write(*,*)'please choose model=2 to have access to tan_Beta'
            STOP
         elseif(model.eq.2) then    
            if(Q_loop.eq.0) then
               if (lglobalprint) write(*,*)'for Higgs-mixing in eff. theory'
               if (lglobalprint) write(*,*)'please choose higgs_mix=1'
               STOP
            else if(Q_loop.eq.1) then
               if (lglobalprint) write(6,*) 'cp-even Higgs (coupling via top-loop) '
               HloopTop      = .true.
               cplModHtop    = cp_even_mod*m_top**2/tan_Beta
               if (lglobalprint) write(6,*) 'cp-odd Higgs (coupling via top-loop) '
               AloopTop      = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)/tan_Beta
            else if(Q_loop.eq.2) then
               if (lglobalprint) write(6,*) 'cp-even Higgs (coupling via bottom-loop) '
               HloopBottom   = .true.
               cplModHbottom = cp_even_mod*(BMASS_H**2)*tan_Beta
               if (lglobalprint) write(6,*) 'cp-odd Higgs (coupling via bottom-loop) '
               AloopBottom   = .true.
               cplModAbottom = cp_odd_mod*(BMASS_H**2)*tan_Beta
            else if(Q_loop.eq.3) then
               if (lglobalprint) write(6,*) 'cp-even SM Higgs: (coupling via top-loop)'
               if (lglobalprint)  write(6,*) 'and (coupling via bottom-loop) '
               HloopTop      = .true.
               HloopBottom   = .true.
               cplModHtop    = cp_even_mod*(m_top**2)/tan_Beta
               cplModHbottom = cp_even_mod*(BMASS_H**2)*tan_Beta
               if (lglobalprint) write(6,*) 'cp-odd Higgs: (coupling via top-loop) '
               if (lglobalprint)  write(6,*) 'and (coupling via bottom-loop) '
               AloopTop      = .true.
               AloopBottom   = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)/tan_Beta
               cplModAbottom = cp_odd_mod*(BMASS_H**2)*tan_Beta
            endif               ! end of Q_loop choice
         endif
      else 
         if (lglobalprint) write(*,*)'Please input HIGGS_MIX = 0,1,2'
         stop
      endif

      if (lglobalprint) write(6,*) '-----------------------------------------------'


     
      END

