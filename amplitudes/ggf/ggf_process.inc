*----------------------------------------------------------------------
c  author   : Michael Kubocz
c  date     : 25.08.2006
c  modified : 24.05.2018 (MK)
*----------------------------------------------------------------------
c modified: Francisco Campanario
c date : 5.11.2012
c add: Variables for Gram determinants and QUAD precision staff
c
ccccccccccccccccccccccccccccccccccccccc
************************** subprocess choice **************************
      logical subprQQ,subprQG,subprGG

******************* Higgs coupling type parameters ********************
      integer Q_loop
      common/qloopg/q_loop
      logical AloopTop,AloopBottom,HloopTop,HloopBottom,
     &        Aeff,Heff,
     &        sHloopTop,sHloopBottom,bHloopTop,bHloopBottom,
     &        sHloopStop1,sHloopSbottom1,bHloopStop1,bHloopSbottom1,
     &        sHloopStop2,sHloopSbottom2,bHloopStop2,bHloopSbottom2
 
*************************** Higgs couplings ***************************
      complex*16 cplModHeff,cplModAeff
      real*8 cplModAtop,cplModAbottom,cplModHtop,cplModHbottom,
     &       cplModsHtop,cplModsHbottom,cplModbHtop,cplModbHbottom,
     &       cp_even_mod,cp_odd_mod,
     &       cplModsHstop(2),cplModsHsbottom(2),
     &       cplModbHstop(2),cplModbHsbottom(2)
      integer higgs_mix
      logical eff_th_corr,eff_th_corr_taylor

*************************** mass parameters ***************************
      real*8 m_top,m_bottom,m_stop(2),m_sbottom(2)

**************************** mixing anlges ****************************
      real*8 sin_Alpha,cos_Alpha,tan_Beta,sin_Beta,cos_Beta

      COMMON/process_settings/cplModAtop,cplModAbottom,cplModHtop,cplModHbottom,
     &      cplModHeff,cplModAeff,
     &      cplModsHtop,cplModsHbottom,cplModbHtop,cplModbHbottom,
     &      cp_even_mod,cp_odd_mod,
     &      m_top,m_bottom,
     &      sin_Alpha,cos_Alpha,tan_Beta,sin_Beta,cos_Beta,
     &      subprQQ,subprQG,subprGG,AloopTop,AloopBottom,HloopTop,HloopBottom,
     &      sHloopStop1,sHloopSbottom1,bHloopStop1,bHloopSbottom1,
     &      sHloopStop2,sHloopSbottom2,bHloopStop2,bHloopSbottom2,
     &      Aeff,Heff,
     &      sHloopTop,sHloopBottom,bHloopTop,bHloopBottom


c      COMMON/process_settings/processGGF,
c     &      subprQQ,subprQG,subprGG,
c     &      AloopTop,AloopBottom,HloopTop,HloopBottom,
c     &      Aeff,Heff,
c     &      sHloopTop,sHloopBottom,bHloopTop,bHloopBottom,
c     &      cplModAtop,cplModAbottom,cplModHtop,cplModHbottom,
c     &      cplModHeff,cplModAeff,
c     &      cplModsHtop,cplModsHbottom,cplModbHtop,cplModbHbottom,
c     &      cp_even_mod,cp_odd_mod,
c     &      m_top,m_bottom,
c     &      mA,
c     &      sin_Alpha,cos_Alpha,tan_Beta,sin_Beta,cos_Beta

c      COMMON/process_settings/subprQQ,subprQG,subprGG,	
c     &      AloopTop,AloopBottom,HloopTop,HloopBottom,
c     &      Aeff,Heff,
c     &      sHloopTop,sHloopBottom,bHloopTop,bHloopBottom,
c     &      cplModAtop,cplModAbottom,cplModHtop,cplModHbottom,
c     &      cplModHeff,cplModAeff,
c     &      cplModsHtop,cplModsHbottom,cplModbHtop,cplModbHbottom,
c     &      cp_even_mod,cp_odd_mod,
c     &      m_top,m_bottom,
c     &      sin_Alpha,cos_Alpha,tan_Beta,sin_Beta,cos_Beta


        integer hexinst(10,2,2),peninst(10,2,2)
        integer boxinst(10,2,2),verinst(10,2,2)
        common/process_insta/hexinst,peninst,boxinst,verinst

        logical HexSingular(2),pensingular(2),boxSingular(2),verSingular(2)
        common/hexSin/HexSingular,pensingular,boxSingular,verSingular

        integer counter(10),gf_iteration,counterC(10),counterD(10)
	logical singular
        COMMON /DetCount/singular,counter,gf_iteration,counterC,counterD

	logical hexeva,peneva,boxeva
	common/evaswitch/hexeva,peneva,boxeva

	logical hexdbleSin, hexQUADSin
	logical pendbleSin, penQUADSin
	logical boxdbleSin, boxQUADSin
	common/pensing/pendbleSin, penQUADSin
	common/boxsing/boxdbleSin, boxQUADSin
	logical verdbleSin, verQUADSin
	common/versing/verdbleSin, verQUADSin



	logical mguu,mguub,mggu,equal_f
	common/mgflavour/mguu,mguub,mggu,equal_f

	integer sub_QG_part,sub_QQ_part
	common/subproces_part/sub_QG_part,sub_QQ_part
