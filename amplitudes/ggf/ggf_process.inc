*----------------------------------------------------------------------
c  author   : Michael Kubocz
c  date     : 25.08.2006
c  modified : 03.03.2011 (MK)
*----------------------------------------------------------------------

************************** subprocess choice **************************
      logical subprQQ,subprQG,subprGG

******************* Higgs coupling type parameters ********************
      integer Q_loop
      logical AloopTop,AloopBottom,HloopTop,HloopBottom,
     &        Aeff,Heff,
     &        sHloopTop,sHloopBottom,bHloopTop,bHloopBottom
 
*************************** Higgs couplings ***************************
      real*8 cplModAtop,cplModAbottom,cplModHtop,cplModHbottom,
     &       cplModHeff,cplModAeff,
     &       cplModsHtop,cplModsHbottom,cplModbHtop,cplModbHbottom,
     &       cp_even_mod,cp_odd_mod
      integer higgs_mix

*************************** mass parameters ***************************
      real*8 m_top,m_bottom

**************************** mixing anlges ****************************
      real*8 sin_Alpha,cos_Alpha,tan_Beta,sin_Beta,cos_Beta

      COMMON/process_settings/cplModAtop,cplModAbottom,cplModHtop,cplModHbottom,
     &      cplModHeff,cplModAeff,
     &      cplModsHtop,cplModsHbottom,cplModbHtop,cplModbHbottom,
     &      cp_even_mod,cp_odd_mod,
     &      m_top,m_bottom,
     &      sin_Alpha,cos_Alpha,tan_Beta,sin_Beta,cos_Beta,
     &      subprQQ,subprQG,subprGG,AloopTop,AloopBottom,HloopTop,HloopBottom,
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


