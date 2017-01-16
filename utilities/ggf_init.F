***********************************************************************
      SUBROUTINE init_gluonfusion
          use readinput
*----------------------------------------------------------------------
c... author   : Michael Kubocz
c... modified : 03.03.2011
***********************************************************************

      IMPLICIT NONE

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
      CALL read_Logical("SUBPRQQ",subprQQ,.true.)
      CALL read_Logical("SUBPRQG",subprQG,.true.)
      CALL read_Logical("SUBPRGG",subprGG,.true.)
c... Mixing between different Higgs bosons
      CALL read_int("HIGGS_MIX",higgs_mix,0)
      if (higgs_mix .eq. 1) then
         CALL read_real("CP_EVEN_MOD",cp_even_mod,1d0)
         CALL read_real("CP_ODD_MOD",cp_odd_mod,1d0)

         select CASE (procID)
         CASE(GFHjj, GFHjj_WW, GFHjj_ZZ_ll,GFHjj_ZZ_lnu)
         CASE DEFAULT
            write(*,*)'  '
            write(*,*)'Sorry, we can not yet include higgs decays'
            write(*,*)'into fermions or photons from gluon fusion with '
            write(*,*)'user-defined mixing between CP-odd and '
            write(*,*)'CP-even Higgs bosons'
            stop
         END select
      end if

      call printnfl(.true.)

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
***********************************************************************

      write(6,*)'--------- PROCESS PARAMETERS IN GLUON FUSION ---------'
      write(6,*)''
      if (abs(xmh-2d0*m_top) .le. 1d0) then
        write(6,*)' -->                  ATTENTION                  <--'
        write(6,*)' -->     threshold effects at m_H = 2*m_top      <--'
        write(6,*)' -->    please use a Higgs mass that fulfills    <--'
        write(6,*)' -->          | m_H - 2*m_top | > 1 GeV          <--'
        write(6,*)' -->                  ATTENTION                  <--'
        write(6,*)''
      endif
      write(6,*)'------- FURTHER PROCESS INFO  -------'
      if(Q_loop.eq.1) write(*,"(T4,A,T22,A,F9.4)")"m_top","= ",m_top
      if(Q_loop.eq.2) then 
         write(*,"(T4,A,T22,A,F9.4)")"m_bottom(m_Higgs)","= ",m_bottom
      endif
      if(Q_loop.eq.3) then 
         write(*,"(T4,A,T22,A,F9.4)")"m_top","= ",m_top
         write(*,"(T4,A,T22,A,F9.4)")"m_bottom(m_Higgs)","= ",m_bottom
      endif
      print*,''
      write(6,*)'Gram det in Cij ->', GRAMCTO0
      write(6,*)'Gram det in Dij ->', GRAMDTO0
      print*,''

      write(6,*) '------------------- PROCESS(ES) ---------------------'

      if (subprQQ) write(6,*) 'QUARK - QUARK SCATTERING '
      if (subprQG) write(6,*) 'QUARK - GLUON SCATTERING ' 
      if (subprGG) write(6,*) 'GLUON - GLUON SCATTERING '
***********************************************************************

      if (higgs_mix .eq. 0) then   ! No higgs mixing: default case

         if (model.eq.1) then   ! SM Higgs
            if(Q_loop.eq.0) then
               print*,''
               write(6,*)'cp-even SM Higgs in the mtop->infinity limit'
               Heff          = .true.
               cplModHeff    = 1d0/3d0 
            else if(Q_loop.eq.1) then
               print*,''
               write(6,*) 'cp-even SM Higgs (coupling via top-loop) '
               HloopTop      = .true.
               cplModHtop    = m_top**2
            else if(Q_loop.eq.2) then
               print*,''
               write(6,*) 'cp-even SM Higgs (coupling via bottom-loop) '
               HloopBottom   = .true.
               cplModHbottom = m_bottom**2
            else if(Q_loop.eq.3) then
               print*,''
               write(6,*) 'cp-even SM Higgs (coupling via top-loop) '
               write(6,*) 'cp-even SM Higgs (coupling via bottom-loop) '
               HloopTop      = .true.
               HloopBottom   = .true.
               cplModHtop    = m_top**2
               cplModHbottom = m_bottom**2
            endif
         endif

         if(model.eq.2) then  ! MSSM, no higgs mixing
*       if(HiggsType.eq.0) then
*         write(*,*)'When working in the MSSM, please choose'
*         write(*,*)'HTYPE = 1,2,3'
*         STOP
*       endif
            if(HiggsType.eq.3) then  ! CP odd higgs in MSSM, no mixing
               if(Q_loop.eq.0) then
                  print*,' '
                  write(6,*)'cp-odd MSSM Higgs A'
                  write(6,*)'in the mtop->infinity limit'
                  Aeff          = .true.
                  cplModAeff    = 0.5d0  
               elseif(Q_loop.eq.1) then
                  print*,' '
                  write(6,*) 'cp-odd MSSM Higgs A:'
                  write(6,*) '(coupling via top-loop) '
                  AloopTop      = .true.
                  cplModAtop    = (m_top**2)/tan_Beta
               elseif(Q_loop.eq.2) then
                  print*,''
                  write(6,*) 'cp-odd MSSM Higgs A:'
                  write(6,*) '(coupling via bottom-loop) '
                  AloopBottom   = .true.
                  cplModAbottom = (m_bottom**2)*tan_Beta
               elseif(Q_loop.eq.3) then
                  print*,''
                  write(6,*) 'cp-odd MSSM Higgs A:'
                  write(6,*) '(coupling via top-loop) '
                  write(6,*) '(coupling via bottom-loop) '
                  AloopTop      = .true.
                  AloopBottom   = .true.
                  cplModAtop    = (m_top**2)/tan_Beta
                  cplModAbottom = (m_bottom**2)*tan_Beta
               endif
            endif

            if(HiggsType.eq.1) then  ! MSSM, light cp-even h0, no higgs mixing
               if(Q_loop.eq.0) then
                  write(*,*)'process is not available'
                  STOP
               endif
               if(Q_loop.eq.1) then
                  print*,''
                  write(6,*) 'cp-even MSSM Higgs h:'
                  write(6,*) '(coupling via top-loop) '
                  sHloopTop     = .true.
                  cplModsHtop   =  -(m_top**2)*cos_alpha/sin_Beta
               elseif(Q_loop.eq.2) then
                  print*,''
                  write(6,*) 'cp-even MSSM Higgs h:'
                  write(6,*) '(coupling via bottom-loop) '
                  sHloopBottom  = .true.
                  cplModsHbottom= (m_bottom**2)*sin_alpha/cos_Beta
               elseif(Q_loop.eq.3) then
                  print*,''
                  write(6,*) 'cp-even MSSM Higgs h:'
                  write(6,*) '(coupling via top-loop) '
                  write(6,*) '(coupling via bottom-loop) '
                  sHloopTop     = .true.
                  sHloopBottom  = .true.
                  cplModsHtop   = -(m_top**2)*cos_alpha/sin_Beta
                  cplModsHbottom= (m_bottom**2)*sin_alpha/cos_Beta 
               endif
            endif

            if(HiggsType.eq.2) then  !  MSSM, cp-even H, no higgs mixing
               if(Q_loop.eq.0) then
                  write(*,*)'process is not available'
                  STOP
               endif
               if(Q_loop.eq.1) then
                  print*,''
                  write(6,*) 'cp-even MSSM Higgs H:'
                  write(6,*) '(coupling via top-loop) '
                  bHloopTop     = .true.
                  cplModbHtop   = -(m_top**2)*sin_alpha/sin_Beta 
               elseif(Q_loop.eq.2) then
                  print*,''
                  write(6,*) 'cp-even MSSM Higgs H:'
                  write(6,*) '(coupling via bottom-loop) '
                  bHloopBottom  = .true.
                  cplModbHbottom= -(m_bottom**2)*cos_alpha/cos_Beta 
               elseif(Q_loop.eq.3) then
                  print*,''
                  write(6,*) 'cp-even MSSM Higgs H:'
                  write(6,*) '(coupling via top-loop) '
                  write(6,*) '(coupling via bottom-loop) '
                  bHloopTop     = .true.
                  bHloopBottom  = .true.
                  cplModbHtop   = -(m_top**2)*sin_alpha/sin_Beta
                  cplModbHbottom= -(m_bottom**2)*cos_alpha/cos_Beta 
               endif
            endif
         endif
      

      else if (higgs_mix .eq. 1) then ! Mixing cp-odd and -even, user defined


         write(6,*)'CP-violating model:'
         write(6,*)'Higgs coupling->(1+i*gamma_5)'
         write(*,"(T4,A,T22,A,F9.4)")"cp_even_mod","= ",cp_even_mod
         write(*,"(T4,A,T22,A,F9.4)")"cp_odd_mod","= ",cp_odd_mod
         print*,''
         
         if(model.eq.1) then    ! SM: higgs = cp_even*H_SM + cp_odd*A_SM
!                                 note: A_SM is defined as A_MSSM with TB=1
            if(Q_loop.eq.0) then
               write(6,*) 'cp-even Higgs in the mtop->infinity limit'
               Heff          = .true.
               cplModHeff    = cp_even_mod*1d0/3d0 
               write(6,*) 'cp-odd Higgs in the mtop->infinity limit'
               Aeff          = .true.
               cplModAeff    = cp_odd_mod*0.5d0  
            else if(Q_loop.eq.1) then
               write(6,*) 'cp-even Higgs (coupling via top-loop) '
               HloopTop      = .true.
               cplModHtop    = cp_even_mod*m_top**2
               write(6,*) 'cp-odd Higgs (coupling via top-loop) '
               AloopTop      = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)
            else if(Q_loop.eq.2) then
               write(6,*) 'cp-even Higgs (coupling via bottom-loop) '
               HloopBottom   = .true.
               cplModHbottom = cp_even_mod*m_bottom**2
               write(6,*) 'cp-odd Higgs (coupling via bottom-loop) '
               AloopBottom   = .true.
               cplModAbottom = cp_odd_mod*(m_bottom**2)
            else if(Q_loop.eq.3) then
               write(6,*) 'cp-even SM Higgs: (coupling via top-loop)'
               write(6,*) 'and (coupling via bottom-loop) '
               HloopTop      = .true.
               HloopBottom   = .true.
               cplModHtop    = cp_even_mod*(m_top**2)
               cplModHbottom = cp_even_mod*(m_bottom**2)
               write(6,*) 'cp-odd Higgs: (coupling via top-loop) '
               write(6,*) 'and (coupling via bottom-loop) '
               AloopTop      = .true.
               AloopBottom   = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)
               cplModAbottom = cp_odd_mod*(m_bottom**2)
            endif               ! end of SM Q_loop choice
            
         else !  MSSM: user-defined mixing, higgs = cp_even*h0 + cp_odd*A0

            if(Q_loop.eq.1) then
               write(6,*) 'cp-odd MSSM Higgs A (via top-loop) '
               AloopTop = .true.
               cplModAtop = cp_odd_mod*(m_top**2)/tan_Beta
               write(6,*) 'cp-even MSSM Higgs h (via top-loop) '
               sHloopTop = .true.
               cplModsHtop = -cp_even_mod*(m_top**2)*cos_alpha/
     -              sin_Beta
            else if(Q_loop.eq.2) then
               write(6,*) 'cp-odd MSSM Higgs A (via bottom-loop) '
               AloopBottom   = .true.
               cplModAbottom = cp_odd_mod*(m_bottom**2)*tan_Beta
               write(6,*) 'cp-even MSSM Higgs h (via bottom-loop) '
               sHloopBottom   = .true.
               cplModsHbottom = cp_even_mod*(m_bottom**2)*sin_alpha/
     -              cos_Beta
            else if(Q_loop.eq.3) then
               write(6,*) 'cp-odd MSSM Higgs A (via top-loop) '
               write(6,*) 'cp-odd MSSM Higgs A (via bottom-loop) '
               AloopTop      = .true.
               AloopBottom   = .true.
               cplModAtop    = cp_odd_mod*(m_top**2)/tan_Beta
               cplModAbottom = cp_odd_mod*(m_bottom**2)*tan_Beta
               write(6,*) 'cp-even MSSM Higgs h (via top-loop) '
               write(6,*) 'cp-even MSSM Higgs h (via bottom-loop) '
               sHloopTop      = .true.
               sHloopBottom   = .true.
               cplModsHtop    = -cp_even_mod*(m_top**2)*cos_alpha/
     -              sin_Beta
               cplModsHbottom = cp_even_mod*(m_bottom**2)*sin_alpha/
     -              cos_Beta
            end if              ! end of Q_loop choice for MSSM
         endif                  ! end of SM/MSSM for user-defined mixing


      else if (higgs_mix .eq. 2) then !  Mixing between MSSM Higgs bosons 
!                                        according to ZHiggs

         if (model .eq. 1) then ! SM
            write(*,*)'Sorry, we can not mix according to the'
            write(*,*)'Z matrix while working in the SM'
            stop
         else                   ! MSSM
            write(*,*)'Higgs mixing according to Z propagator matrix:'

            if (Q_loop .eq. 0) then
               write(*,*)'Sorry, we can not mix Higgs bosons '
               write(*,*)'according to the Z matric with Q_loop = 0'
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
               cplModsHbottom= DBLE(ZHiggs(HiggsType,1))*(m_bottom**2)*
     -              sin_alpha/cos_Beta
               bHloopBottom  = .true.
               cplModbHbottom= -DBLE(ZHiggs(HiggsType,2))*(m_bottom**2)*
     -              cos_alpha/cos_Beta 
               AloopBottom   = .true.
               cplModAbottom = DBLE(ZHiggs(HiggsType,3))*(m_bottom**2)*
     -              tan_Beta
            else if(Q_loop.eq.3) then
               sHloopTop     = .true.
               cplModsHtop   = -DBLE(ZHiggs(HiggsType,1))*(m_top**2)*
     -              cos_alpha/sin_Beta
               sHloopBottom  = .true.
               cplModsHbottom= DBLE(ZHiggs(HiggsType,1))*(m_bottom**2)*
     -              sin_alpha/cos_Beta 
               bHloopTop     = .true.
               cplModbHtop   = -DBLE(ZHiggs(HiggsType,2))*(m_top**2)*
     -              sin_alpha/sin_Beta
               bHloopBottom  = .true.
               cplModbHbottom= -DBLE(ZHiggs(HiggsType,2))*(m_bottom**2)*
     -              cos_alpha/cos_Beta 
               AloopTop      = .true.
               cplModAtop    = DBLE(ZHiggs(HiggsType,3))*(m_top**2)/
     -              tan_Beta
               AloopBottom   = .true.
               cplModAbottom = DBLE(ZHiggs(HiggsType,3))*(m_bottom**2)*
     -              tan_Beta
            end if              ! end of Q_loop choice
            
         end if                 ! end of SM/MSSM loop
         
      else 
         write(*,*)'Please input HIGGS_MIX = 0,1,2'
         stop
      endif


      write(6,*) '-----------------------------------------------------'

     
      END

